# cython: language_level=3str, boundscheck=False, wraparound=False

from warnings import warn

from numpy import absolute, amax, arange, array, cos, copy, empty, exp, eye, zeros, sqrt as np_sqrt, savez, ones, outer, sin as np_sin, sum, mean, pi as np_pi, NaN, Inf
from scipy.linalg import det, svd
from scipy.optimize import minimize
from libc.math cimport sin, sinh, sqrt, pi, fabs, asin
from numpy.linalg import solve, LinAlgError
from scipy.optimize import minimize
from time import sleep
from tools import absq
from pickle import dump as pkl_dump, load as pkl_load, HIGHEST_PROTOCOL

I5I2tme = 2 * (pi/2) ** 3 / (3 * sqrt(15))

class Grapher:
    def __init__(self, samplingSteps = 128):
        self.sL = arange(samplingSteps) / samplingSteps

        self.alphaX, self.betaX, self.gammaX, self.betaY, self.k, \
            self.eta, self.etaPrime, self.b, self.bSq, self.bCub, = [empty(samplingSteps) for _ in range(10)]

        self.mSext = zeros(samplingSteps)
        self._B = empty((2, samplingSteps))

    def feedX(self, homX, homXPrime):
        sqIx = np_sqrt((homX.conj() * homXPrime).imag) # this 'invariant' is technically not invariant in s
        homX /= sqIx
        homXPrime /= sqIx
        self.betaX[:] = absq(homX)
        self.gammaX[:] = absq(homXPrime)
        self.alphaX[:] = -(homX.conj() * homXPrime).real

    def feedY(self, homY, homYPrime):
        self.betaY[:] = absq(homY)
        self.betaY /= absolute((homY.conj() * homYPrime).imag)

    def feedB(self, b):
        self.b = b
        self.bSq = b ** 2
        self.bCub = self.bSq * absolute(self.b)

    def i2(self):
        return sum(self.bSq)

    def i3(self):
        return sum(self.bCub)

    def i4(self):
        # these are sums not means
        return 2*sum(self.eta * self.b * self.k)

    def i5(self):
        return sum(self.bCub * (self.gammaX * self.eta**2 + 2 * self.alphaX * self.eta * self.etaPrime
                                + self.betaX * self.etaPrime ** 2))

    def F(self):
        i2 = self.i2()
        i4 = self.i4()
        # if absolute(i4) > i2:
        #     return Inf
        return self.i5() / (i2 - i4) / I5I2tme

    def jX(self):
        return 1 - self.i4() / self.i2()

    def naturalChroma(self):
        chroma = empty(2)
        chroma[0] = -mean(self.k * self.betaX)  # x plane
        chroma[1] =  mean(self.k * self.betaY)  # y plane
        chroma /= 4
        return chroma

    def fullChroma(self):
        mEtaK = self.k - self.mSext * self.eta
        chroma = empty(2)
        chroma[0] = -mean(self.betaX * mEtaK)
        chroma[1] =  mean(self.betaY * mEtaK)
        chroma /= 4
        return chroma

    def _matrixB(self):
        self._B[0] = -self.betaX * self.eta
        self._B[1] =  self.betaY * self.eta
        self._B /= 4*self.sL.size
        return self._B


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

cdef _fillB(double[:] b, double[:] b_nonzero, int P):
    cdef int pbmax = b_nonzero.size+1
    cdef int p
    for p in range(1,pbmax):
        b[P-p] = b[p+P] = b_nonzero[p-1]


class FloquetSolve:
    def __init__(self, P=32):
        self.P = P
        self.pRange = arange(-P, P + 1)
        self._mat0 = eye(self.pRange.size)
        self._mat = eye(self.pRange.size)
        self._denom0 = empty(2*P+1)
        self.b = zeros(2*P+1)
        self.b[P] = 1
        self.k = None

    # -- setting k and obtaining linear optics --
    # as .setK can be used for multiple functions .C, .tune and .tuneUPars, calling it returns no result

    def setK(self, k):
        self.k = k
        self._denom0 = self.k[0] - 4 * (self.pRange) ** 2
        vals_np = self.k[:,None] / self._denom0[None,:]
        populateM(self._mat0, vals_np)

    def C(self):
        # Whittaker-Hill formula for tune estimates,
        # eq. (2.13) in refs/characteristicExponents.pdf
        return _C(self.k[0], det(self._mat0))

    def tune(self):
        return _tune(self.k[0], det(self._mat0))

    def tuneUPars(self):
        tune = _tune(self.k[0], det(self._mat0))
        denom = self.k[0] - 4 * (self.pRange + tune) ** 2
        vals_np = self.k[:,None] / denom[None,:]
        populateM(self._mat, vals_np)
        try:
            u, s, vh = svd(self._mat)
            uPars = vh[-1, :].T  # from null_space
            err = sum(absq(self._mat @ uPars))
            if err > 1e-10:
                warn(': high error at tuneX=%.2f' % tune)
        except ValueError:
            raise ValueError('SVD collapsed at tuneX=%.2f' % tune)
        return tune, uPars

    # -- end of setting k --

    def setB(self, bNonZero):
        _fillB(self.b, bNonZero, self.P)
        vPars = solve(self._mat0, self.b / self._denom0)
        return vPars


class FourierCell:
    """contains functions that connect .gr and .fs or that are helper functions."""
    def __init__(self):
        self.k = ones(2, float) / 2
        self.bNonZero = ones(2, float) / 2

        self.gr = Grapher()
        self.fs = FloquetSolve()

        self._convert = exp(2.j*np_pi*outer(self.gr.sL, self.fs.pRange))
        # this must be update when k.size changes:
        self._convertK = copy(self._convert[:,self.fs.P:(self.fs.P+self.k.size)].real)
        self._convertK[:,1:] *= 2  # k_0 + 2 k_1 cos + 2 k_2 cos...

        self.desiredTune = array([-1,-1], dtype=float)
        self._desiredC = array([-1,-1], dtype=float)

    def setKx(self, k):
        # compute tune, alpha, beta, gamma
        self.fs.setK(k)
        tune, uParsX = self.fs.tuneUPars()

        #euSpiral = exp(2.j * tune * np_pi * self.gr.sL)
        uX = (self._convert @ uParsX) #* euSpiral
        uXPrime = 2.j * (self._convert @ (uParsX * (self.fs.pRange + tune))) #* euSpiral
        self.gr.feedX(uX, uXPrime)

        self.gr.k[:] = self._convertK @ k
        self.k = k

    def setKxy(self, k):
        self.fs.setK(-k)
        tuneY, uParsY = self.fs.tuneUPars()

        #euSpiral = exp(2.j * tuneY * np_pi * self.gr.sL)
        uY = (self._convert @ uParsY) #* euSpiral
        uYPrime = 2.j * self._convert @ (uParsY * (self.fs.pRange + tuneY)) #* euSpiral
        self.gr.feedY(uY, uYPrime)
        self.setKx(k)

    def setB(self, bNonZero):
        """use after .setKx or .setKxy"""
        vPars = self.fs.setB(bNonZero)
        self.gr.feedB(self._convert.real @ self.fs.b)
        self.gr.eta[:] = self._convert.real @ vPars
        self.gr.etaPrime[:] = -2 * self._convert.imag @ (vPars * self.fs.pRange)
        self.bNonZero = bNonZero

    def lagrangeTuneEq(self, k):
        self.fs.setK(-k)
        f = (self.fs.C() - self._desiredC[1]) ** 2
        self.fs.setK(k)
        f += (self.fs.C() - self._desiredC[0]) ** 2
        return f

    def tuneTo(self, tuneX, tuneY):
        self.desiredTune[:] = (tuneX, tuneY)
        self._desiredC = np_sin(np_pi * self.desiredTune) ** 2

        kStart=copy(self.k)

        res = minimize(self.lagrangeTuneEq, kStart)
        self.setKxy(res.x)

    def mWedge2(self, setMSext=False):
        # compensate chroma by sextupole vector with 2 dof. only implemented for k.size=2
        twoTimesTwo = self.gr._matrixB() @ self._convertK
        try:
            x = solve(twoTimesTwo,self.gr.naturalChroma())
        except LinAlgError:
            raise LinAlgError('k.size != 2?')
        if setMSext:
            self.gr.mSext[:] = self._convertK @ x
        x = absolute(x)
        return x[0] + 2*x[1]

    def G2(self, setMSext=False):
        return self.gr.F() * self.mWedge2(setMSext)**(0.75)

    def mWedgeN(self, basisConstraints, chroma):
        # minimization of sextupole excursion with arbitrary dof by linear programming
        # linprog(...)
        # return 1/sol[0]
        pass

    def __str__(self):
        return "FourierCell at k=%s. graphed F=%.2f, jX=%.2f" % (self.k.__str__(), self.gr.F(), self.gr.jX())

def directSearch(fc, fun):
    b1=array([0.0])
    step=-0.2 #-0.2
    oldVal=Inf
    while absolute(step)>0.002: #0.005
        b1+=step
        if b1<-5:
            break
        fc.setB(b1)
        newVal=fun()
        if newVal>oldVal and oldVal>0:
            step /= -2
        oldVal=newVal
    return b1, oldVal


class Map:
    def __init__(self, shape, name='', atNames=None):
        self.name = name
        self.b1 = empty(shape)
        self.fval = empty(shape)
        self.atNames = atNames
        self.atArray = empty((*shape,len(atNames))) if atNames is not None else None

    def write(self, q, p, b1, fval):
        self.b1[q, p] = b1
        self.fval[q, p] = fval

    def save(filename):
        with open(filename,'w') as f:
        pkl_dump(self, f, HIGHEST_PROTOCOL)

    def load(filename):
        with open(filename,'r') as f:
        self = pkl_load(f)


class TuneMap:
    def __init__(self, dots=100, tuneRange=None):
        if tuneRange is None:
            tuneRange = 0.5 * (arange(1, dots) / dots)
        self.tuneX = tuneRange[None, :]
        self.tuneY = tuneRange[:, None]

        shape = (self.tuneY.size, self.tuneX.size)
        self.k = empty((*shape,2))
        self.k[:] = NaN
        self.chroma = copy(self.k)

        self.mapF = Map(shape, name='F', atNames=('jX',))
        self.mapG2 = Map(shape, name='G2', atNames=('jX', 'F', 'mWedge2'))

    def make(self):
        fc = FourierCell()
        for q, tuneY in enumerate(self.tuneY.flat):
            print(q)
            for p, tuneX in enumerate(self.tuneX.flat):
                fc.tuneTo(tuneX, tuneY)
                self.k[q,p] = fc.k
                self.chroma[q,p] = fc.gr.naturalChroma()
                # print('tX=%.2f, ty=%.2f' % (tuneX, tuneY))
                # ToDo:
                # - write newton search on Jx (start, end) first
                self.mapF.write(q,p,*directSearch(fc, fc.gr.F))
                self.mapF.atArray[q,p,0] = fc.gr.jX()
                # print('G2')
                self.mapG2.write(q,p,*directSearch(fc, fc.G2))
                self.mapG2.atArray[q,p] = fc.gr.jX(), fc.gr.F(), fc.mWedge2(False)
        self.mapF.save('F.pkl')
        self.mapG2.save('G2.pkl')
        savez('tunechroma.npz', tuneX=self.tuneX, tuneY=self.tuneY, k=self.k, chroma=self.chroma)
