from numpy import mean, zeros, arange, pi, absolute, sin, sqrt, empty, \
        exp, NaN, empty_like, moveaxis, amin, array, \
        amax, polyfit, polyval, cos, ones_like, zeros_like, ones, copy
from numpy.linalg import LinAlgError
from numpy.linalg import solve
from scipy.optimize import minimize
from tools import loadData, saveData, centaur

from microMacroCell import SimpleFloquetCell, OptFloquetCell, Grapher, I5I2tme


class FloquetKSpace:
    def __init__(self, *kRanges):
        shape = [kRange.size for kRange in kRanges]
        self.kRanges = kRanges

        b1range=zeros(1, float)
        self.whm = OptFloquetCell(kDim=len(shape))
        self.cX = empty(shape, float)
        self._cY = empty_like(self.cX)
        self.withY = False

    def solveCx(self):
        k = empty(2, float)
        for p, k[0] in enumerate(self.kRanges[0]):
            for q, k[1] in enumerate(self.kRanges[1]):
                self.cX[p, q] = self.whm.C(k)

    def solveCxy(self):
        k = empty(2, float)
        for p, k[0] in enumerate(self.kRanges[0]):
            for q, k[1] in enumerate(self.kRanges[1]):
                self.cX[p, q] = self.whm.C(k)
                self._cY[p, q] = self.whm.C(-k)
        self.withY = True

    @property
    def cY(self):
        if self.withY:
            return self._cY
        elif (self.kRanges[0][0] + self.kRanges[0][-1]) > 1e-10:
            raise Exception('symmetric k0range required as solveCxy not used')
        else:
            return self.cX[::-1]

    def plotStability(self, ax, axRanges=None, withMirror=True, labels=('$k_0$','$k_1$'), borderLine=0.7, borderColors=('blue', 'green'), tuneLevels=None):
        boundaries = 0, 1
        levelColor = '0.4'
        if axRanges is None:
            axRanges = self.kRanges[:2]
        ax.contourf(axRanges[0], axRanges[1], self.cX.T, boundaries, colors='#ccccff')
        ax.contour(axRanges[0], axRanges[1], self.cX.T, boundaries, colors=('black', borderColors[0]), linewidths=borderLine)
        if tuneLevels is not None:
            levels = sin(pi * tuneLevels) ** 2  # before [0, 0.25, 0.5, 0.75, 1]
            ax.contour(axRanges[0], axRanges[1], self.cX.T, levels, colors=levelColor, linewidths=0.5)
        if withMirror:
            cY = self.cY
            ax.contourf(axRanges[0], axRanges[1], cY.T, boundaries, colors='#ccffcc88')
            ax.contour(axRanges[0], axRanges[1], cY.T, boundaries, colors=('black', borderColors[1]), linewidths=borderLine)
            if tuneLevels is not None:
                ax.contour(axRanges[0], axRanges[1], cY.T, levels, colors=levelColor, linewidths=0.5)
        ax.set(xlabel=labels[0], ylabel=labels[1])


class TuneMap:
    def __init__(self, tuneRange=None, dots=100):
        if tuneRange is None:
            tuneRange = 0.5 * (arange(1, dots) / dots)
        self.tuneX = tuneRange[None, :]
        self.tuneY = tuneRange[:, None]

        shape = (self.tuneY.size, self.tuneX.size)
        self.k0 = empty(shape)
        self.k1 = empty(shape)
        self.b1 = empty(shape)
        self.k0[:] = self.k1[:] = self.b1[:] = NaN

        self.cell = OptFloquetCell()

    def generateKs(self):
        k = 0.5*ones(2, float)
        for q, tuneY in enumerate(self.tuneY.flat):
            print(q)
            for p, tuneX in enumerate(self.tuneX.flat):
                k = self.cell.tuneTo(tuneX, tuneY, k)
                self.k0[q, p] = k[0]
                self.k1[q, p] = absolute(k[1])

    def saveKs(self, npz_filename):
        saveData(npz_filename, tuneX=self.tuneX, tuneY=self.tuneY, k0=self.k0, k1=self.k1, b1=self.b1)

    def loadKs(self, npz_filename):
        x = loadData(npz_filename)
        self.tuneX, self.tuneY, self.k0, self.k1, self.b1 = [x[fld] for fld in ('tuneX', 'tuneY', 'k0', 'k1', 'b1')]

    def mapper(self, whm_method, nargout=1, dtype=float):
        outputs = empty((*self.k0.shape, nargout), dtype)
        def obifun(b1, k0, k1):
            x = array([k0,k1,b1], dtype=float) # br note: b -> b1
            return whm_method(x)

        for q in range(self.tuneY.size):
            print(q)
            for p in range(self.tuneX.size):
                outputs[q, p] = minimize(obifun, 0.0, args=(self.k0[q, p], self.k1[q,p]))
        return moveaxis(outputs, self.k0.ndim, 0)

    def generateCubes(self):
        (tuneX, tuneY, chromaX, chromaY, self.cubes['I4'], self.cubes['I5'], self.cubes['momComp'],
         self.cubes['m0'], self.cubes['m1']) = self.mapper(self.cell.scalarSolve, nargout=9)
        self.cubes.update({'tuneX': tuneX[:, :, None], 'tuneY': tuneY[:, :, None],
                           'chromaX': chromaX[:, :, None], 'chromaY': chromaY[:, :, None]})

    def saveCubes(self, npz_cubefile):
        saveData(npz_cubefile, **self.cubes)

    def loadCubes(self, npz_cubefile):
        self.cubes = dict(loadData(npz_cubefile))

    def processCubes(self):
        F = self.cubes['I5'] / (self.cubes['I2'] - self.cubes['I4'])
        F /= I5I2tme
        Jx = 1 - self.cubes['I4'] / self.cubes['I2']
        F[Jx < 0] = 1e8
        # F[Jx > 2] = 1e8

        mWedge1 = absolute(self.cubes['m0']) + 2 * absolute(self.cubes['m1'])
        G1 = F*mWedge1**(3/4)

        self.cubes.update({'F': F, 'Jx': Jx, 'mWedge1': mWedge1, 'G1': G1})

    def minInB1(self, minCubeName, *otherCubeNames, classic=False):
        minCube = self.cubes[minCubeName]
        otherCubes = [self.cubes[otherCubeName] for otherCubeName in otherCubeNames]
        idx = minCube.argmin(axis=-1)
        optB1 = empty_like(self.k0, float)
        optMinSlice = empty_like(self.k0, float)
        otherSlices = [empty_like(self.k0, float) for _ in range(len(otherCubes))]
        otherRange = arange(len(otherCubes))
        mi, ma = amin(self.b1), amax(self.b1)
        if classic:
            idx = minCube.argmin(axis=-1)
            for q in range(optB1.shape[0]):
                for p in range(optB1.shape[1]):
                    ix = idx[q,p]
                    optB1[q, p] = self.b1[ix]
                    optMinSlice[q, p] = minCube[q,p,ix]
                    for n in otherRange:
                        otherSlices[n][q, p] = otherCubes[n][q, p, ix]
        else: # quadratic interpolation between points
            for q in range(optB1.shape[0]):
                for p in range(optB1.shape[1]):
                    ix = idx[q, p]
                    if ix == 0:
                        continue
                    ply = polyfit(self.b1[ix - 1:ix + 2], minCube[q, p, ix - 1:ix + 2], 2)
                    optB1[q, p] = - ply[1] / (2 * ply[0])
                    optMinSlice[q, p] = ply[2] - (ply[1] ** 2) / (4 * ply[0])
                    if optB1[q, p] < mi or optB1[q, p] > ma:
                        optB1[q, p] = self.b1[-1]
                        optMinSlice[q, p] = minCube[q, p, -1]

                    for n in otherRange:
                        try:
                            ply = polyfit(self.b1[ix - 1:ix + 2], otherCubes[n][q, p, ix - 1:ix + 2], 2)
                            otherSlices[n][q, p] = polyval(ply, optB1[q, p])
                        except LinAlgError:
                            otherSlices[n][q, p] = NaN
        return optB1, optMinSlice, {key: value for (key, value) in zip(otherCubeNames, otherSlices)}

    def grayDiagram(self, ax, quantity, levels=arange(11), fmt='%i', grayDiv=5, **kwargs):
        centaur(ax, self.tuneX.flat, self.tuneY.flat, quantity, levels, fmt=fmt, grayDiv=grayDiv, **kwargs)
        ax.set(xlabel=r'$\nu_x$', ylabel=r'$\nu_y$', xlim=(0, 0.5), ylim=(0, 0.5))
        ax.label_outer()

    def plotJx(self, ax, Jx):
        self.grayDiagram(ax, Jx, levels=arange(0,2.05,0.1), fmt='%.1f', grayDiv=2, faceLims=((2,10),), faceColors=('0.6',))

    def plotF(self, ax, F):
        # plot TME emittance ratio surface
        self.grayDiagram(ax, F, grayMax=2, faceLims=((0,2), (0,1)), faceColors=('#eeffee','#eeeeff'))

    def plotMWedge(self, ax, mWedge):
        self.grayDiagram(ax, mWedge, levels=arange(0,5,0.5), fmt='%.1f', colors='magenta', grayMax=2, grayDiv=5)
