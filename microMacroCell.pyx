# cython: language_level=3str, boundscheck=False, wraparound=False

from warnings import warn

from numpy import absolute, amax, arange, cos, empty, exp, eye, zeros, sqrt as np_sqrt, ones, outer, sum, mean
from scipy.linalg import det, svd
from scipy.optimize import minimize
from libc.math cimport sin, sinh, sqrt, pi, fabs, asin
from twist.tools import absq
from numpy.linalg import solve

cdef void populateM(double[:,:] M, double[:,:] vals):  #, int P):
    cdef int R = M.shape[0]  #=2P+1
    cdef int maxind = M.shape[1]
    cdef int N = vals.shape[0]
    cdef int n, r, q
    for n in range(1,N):
        #for p in range(-P,P+1):
        for r in range(R):
            if r >= n:
                M[r, r-n] = vals[n, r]
            q = r+n
            if q < maxind:
                M[r, q] = vals[n, r]

cdef double _C(double k0, double detM):
    cdef double c = pi * sqrt(fabs(k0)) / 2
    if k0 < 0:
        c = sinh(c)
        c *= -c
    else:
        c = sin(c)
        c *= c
    return c*detM


cdef double _tune(double k0, double detM):
    cdef double c = _C(k0, detM)
    return asin(sqrt(c)) / pi

I5I2tme = 2 * (pi/2) ** 3 / (3 * sqrt(15))

class SimpleFloquetCell:
    # ToDo: now we are on the right track!
    def __init__(self, P=32, kDim=2):
        self.P = P
        self.pRange = arange(-P, P + 1)
        self._mat = eye(self.pRange.size)  # worker array

    def setM(self, k, nu=0.0):
        self._denom = k[0] - 4 * (self.pRange + nu) ** 2
        vals_np = k[:,None] / self._denom[None,:]
        populateM(self._mat, vals_np)

    def C(self, k):
        # Whittaker-Hill formula for tune estimates,
        # eq. (2.13) in refs/characteristicExponents.pdf
        self.setM(k)
        return _C(k[0], det(self._mat))

    def tune(self, k):
        self.setM(k)
        return _tune(k[0], det(self._mat))

    def uPars(self, k, tune):
        self.setM(k, tune)
        try:
            u, s, vh = svd(self._mat)
            uPars = vh[-1, :].T  # from null_space
            err = sum(absq(self._mat @ uPars))
            if err > 1e-10:
                warn(': high error at tuneX=%.2f' % tune)
        except ValueError:
            raise ValueError('SVD collapsed at tuneX=%.2f' % tune)
        return uPars

    def tuneTo(self, tuneX, tuneY, kStart):
        desiredC = [sin(pi * tune) ** 2 for tune in (tuneX, tuneY)]

        def objectiveFunction(k):
            f = (self.C(k) - desiredC[0]) ** 2
            f += (self.C(-k) - desiredC[1]) ** 2
            return f

        res = minimize(objectiveFunction, kStart)
        return res.x


cdef fillB(double[:] b, double[:] b_nonzero, int P):
    cdef int pbmax = b_nonzero.size+1
    cdef int p
    for p in range(1,pbmax):
        b[P-p] = b[p+P] = b_nonzero[p-1]


class Grapher:
    def __init__(self, pRange, int samplingSteps, bSize=1):
        self.s_L = arange(samplingSteps) / samplingSteps
        self.pRange = pRange
        self._eUp = exp(2.j * pi * outer(self.pRange, self.s_L))
        self._cosTwo = 2*cos(2*pi*self.s_L)
        self._cosFour = 2*cos(4*pi*self.s_L)

        self.alphaX, self.betaX, self.gammaX, self.betaY, self.k = [empty(samplingSteps, float) for _ in range(5)]
        self.eta, self.etaPrime, self.b, self.bSq, self.bCub = [empty(samplingSteps, float) for _ in range(5)]

        # self._wedgeMat = empty(?,samplingSteps)
        # self._wedgeMat[2:] = ?

    def inS_lim(self, km):
        # works for k, m (but not b)
        kins = 2 * self._eUp[:, (1 - self.pRange[0]):(km.size - self.pRange[0])].real @ km[1:]
        kins += km[0]
        return kins

    def feedX(self, double tune, uPars, k, vPars):
        eUpp = self._eUp * exp(2.j * tune * pi * self.s_L)
        homX = eUpp @ uPars
        uParPrime = uPars * (self.pRange + tune)
        homXPrime = 2.j * eUpp @ uParPrime
        print(uParPrime.shape)
        sqIx = np_sqrt((homX.conj() * homXPrime).imag) # this 'invariant' is technically not invariant in s
        homX /= sqIx
        homXPrime /= sqIx
        self.betaX[:] = absq(homX)
        self.gammaX[:] = absq(homXPrime)
        self.alphaX[:] = -(homX.conj() * homXPrime).real
        self.k[:] = self.inS_lim(k)
        self.eta[:] = self._eUp.real @ vPars
        self.etaPrime[:] = -2 * self._eUp.imag @ (vPars * self.pRange)

    def feedY(self, double tune, uPars):
        eUpp = self._eUp * exp(2.j * tune * pi * self.s_L)
        homY = eUpp @ uPars
        homYPrime = 2.j * eUpp @ (uPars * (self.pRange + tune))
        self.betaY[:] = absq(homY) / (homY.conj() * homYPrime).imag

    def feedB(self, bPars):
        self.b[:] = self._eUp.real @ bPars
        self.bSq[:] = self.b ** 2
        self.bCub[:] = self.bSq * absolute(self.b)

    def i2(self):
        return sum(self.bSq)

    def i3(self):
        return sum(self.bCub)

    def i45(self):
        # these are sums not means
        return (2*sum(self.eta * self.b * self.k),
               sum(self.bCub * (self.gammaX * self.eta**2 + 2 * self.alphaX * self.eta * self.etaPrime
                                + self.betaX * self.etaPrime ** 2)))

    def momComp(self):
        return mean(self.eta * self.b)

    def chroma(self):
        chromaX = - mean(self.k * self.betaX) / 4
        chromaY = mean(self.k * self.betaY) / 4
        return chromaX, chromaY

    def compensateM01(self, chromaX, chromaY):
        A_00 = mean(self.betaX*self.eta)
        A_10 = mean(self.betaY*self.eta)
        A_01 = mean((self._cosTwo*self.betaX)*self.eta)
        A_11 = mean((self._cosTwo*self.betaY)*self.eta)
        pre = 4 / (A_00 * A_11 - A_10*A_01)
        m0 = -pre*(A_11*chromaX + A_01*chromaY)
        m1 = pre*(A_00*chromaY + A_10*chromaX)
        return m0, m1

    def solveM(self):
        pass

cdef (double, double) _sfjx(double i2, double i4, double i5):
    return i5/(i2-i4), 1 - i4/i2


class OptFloquetCell(SimpleFloquetCell):
    def __init__(self, P=32, kDim=2):
        SimpleFloquetCell.__init__(self, P, kDim)
        self._bPars = zeros(self.pRange.size, float)
        self._bPars[P] = 1
        self.grapher = Grapher(self.pRange, P*16)
        self.kDim = kDim

    def linear(self, k, b_nonzero):
        # vertical plane first, as its not used afterwards
        # some of it follows from symmetry for Mathieu but not in general
        kY = -k
        tuneY = self.tune(kY)
        uParsY = self.uPars(kY, tuneY)

        tuneX = self.tune(k)  #._mat set for nu=0
        fillB(self._bPars, b_nonzero, self.P)
        vPars = solve(self._mat, self._bPars / self._denom) #must be before .uPars

        uParsX = self.uPars(k, tuneX) #._mat set for nu = tuneX, not used further
        return tuneX, tuneY, uParsX, uParsY, vPars

    def scalarSolve(self, k, b_nonzero):
        tuneX, tuneY, uParsX, uParsY, vPars = self.linear( k, b_nonzero )

        # k, b
        self.grapher.feedB(self._bPars)
        self.grapher.feedX(tuneX, uParsX, k, vPars)
        self.grapher.feedY(tuneY, uParsY)
        i2 = self.grapher.i2()
        i4r, i5 = self.grapher.i45()
        # compute chromaticity and compensation
        chromaX, chromaY = self.grapher.chroma()
        m0, m1 = self.grapher.compensateM01(chromaX, chromaY)
        return tuneX, tuneY, i2, i4r, i5, m0, m1 #, sextupole & chromaticity quantities

    # def scaledF(self, x):
    #     tuneX, tuneY, i2Sum, i4rSum, i5Sum, m0, m1 = self.scalarSolve(x[:self.kDim], b_nonzero)
    #     sf, jx = _sfjx(i2Sum,i4rSum,i5Sum)
    #     if jx<0:
    #         sf = 1000.0
    #     print('k0=%.3f,k1=%.3f,b1=%.3f. tX=%.3f,tY=%.3f. F=%.2f,Jx=%.2f.' % (*x, tuneX, tuneY, sf/I5I2tme,jx))
    #     return sf

    def scaledG(self, x):
        # build objective function. input: m_subspace -> map onto (plane & then s space) -> get m_wedge
        tuneX, tuneY, i2Sum, i4rSum, i5Sum, m0, m1 = self.scalarSolve(x[:self.kDim], x[self.kDim:])
        g, jx = _sfjx(i2Sum,i4rSum,i5Sum)
            # solve the m_wedge optimization problem
        mWedge = absolute(m0)+2*absolute(m1)
        g *= mWedge**0.75
        if tuneY < 0.05 or jx<0:
            g = 1000.0
            print('k0=%.3f,k1=%.3f,b1=%.3f. tX=%.3f,tY=%.3f. G=%.2f,Jx=%.2f.' % (*x, tuneX, tuneY, g/I5I2tme,jx))
        return g
