"""Figure generator for Mathieu cells"""
from matplotlib import rc
from matplotlib.pyplot import setp
from numpy import arange, squeeze, array, empty_like, amax, absolute, NaN
from scipy.optimize import minimize as minix, minimize_scalar as minix_scalar
from tools import subplots, saveFig, centaur
from floquetCell import FloquetKSpace
from fourierCell import TuneMap, FourierCell
from time import time
import os

rc('text', usetex=True)
rc('font', family='serif')
# rc('contour', negative_linestyle='solid')

def generateMapIfNotExisting():
    tm = TuneMap()  # b1Range=-arange(0, 1.5, 0.02)) # -arange(0, 1.5, 0.05)  # arange(0, 2, 0.1), arange(-1.4, 0.01, 0.05)
    if os.path.isfile('F.pkl'):
        tm.load()
    else:
        print("computing k values for tune map. this only needs to be done once.")
        tic = time()
        tm.make()
        toc = time()
        print("elapsed time: %.1f sec." % (toc-tic))
    return tm


def grayDiagram(ax, tm: TuneMap, quantity, levels=arange(11), fmt='%i', grayDiv=5, borderColors=('blue', 'green'),
                **kwargs):
    centaur(ax, tm.tuneX.flat, tm.tuneY.flat, quantity, levels, fmt=fmt, grayDiv=grayDiv, **kwargs)
    ax.set(xlabel=r'$\nu_x$', ylabel=r'$\nu_y$', xlim=(0, 0.5), ylim=(0, 0.5))
    ax.spines['right'].set_color(borderColors[0])
    ax.spines['top'].set_color(borderColors[1])
    ax.label_outer()


def plotF(ax, tm, F):
    grayDiagram(ax, tm, F, arange(6), fmt='%i', grayMax=2, grayDiv=10,
                faceLims=((0.0, 1.0),), faceColors=('#ffffcc',))


def plotJx(ax, tm : TuneMap, jX, grayMax=3.0):
    grayDiagram(ax, tm, jX, arange(0, 3.1, 0.5), fmt='%.1f', grayMax=grayMax,
                faceLims=((-10, 0), (3, 100)), faceColors=('#cccccc', '#cccccc'))


def plotMultipoles(ax, fc : FourierCell, sVar='s', sextupoles=False):
    ax.plot(2*fc.gr.sL, fc.gr.b, label='$b(%c)$' % sVar, linewidth=0.7)
    ax.plot(2*fc.gr.sL, fc.gr.k, label='$k(%c)$' % sVar, linewidth=0.7)
    if sextupoles:
        axB.plot(2*fc.gr.sL, fc.gr.mSext, label='$m(%c)$' % sVar, linewidth=0.7)
    ax.axhline(0, color='black', linestyle='dashed', linewidth=0.5)
    ax.set(xlim=(0,1), xlabel=r'%c / $2\pi$' % sVar)


def obig(b2, b1, fc : FourierCell):
    fc.setB(array([b1,b2]))
    return fc.G()  # (fc.gr.jX()-2.5)**2


def b1scan(axA, axB, nuX=0.45, nuY=0.35, minim=False):
    fc = FourierCell()
    tunes = fc.tuneTo(nuX, nuY)
    # tunes = fc.tuneTo(0.15, 0.35)
    print("k_0 = %.5f, k_1 = %.5f" % tuple(fc.k))
    b1range = -arange(0,2.5,0.02)
    F = empty_like(b1range)
    jX = empty_like(b1range)
    i5i2 = empty_like(b1range)
    G = empty_like(b1range)
    maxM = empty_like(b1range)
    for n, b1 in enumerate(b1range):
        if minim:
            # result = minix(obig, 0, args=(b1,fc))
            result = minix_scalar(obig, bounds=(-2,2), args=(b1,fc))
            fc.setB(array([b1,result.x]))
        else:
            fc.setB(array([b1]))
        F[n] = fc.gr.F()
        G[n] = fc.G()
        jX[n] = fc.gr.jX()
        i5i2[n] = fc.gr.i5() / fc.gr.i2()
        maxM[n] = amax(absolute(fc.gr.mSext))
        # print('b1 = %.2f : chrom.nat.x=%.4f,y=%.4f, chrom.full.x=%.4f,y=%.4f' % (b1, *fc.gr.naturalChroma(), *fc.gr.fullChroma()))
    F[jX <= 0] = NaN
    G[jX <= 0] = NaN

    axA.plot(-b1range, i5i2, label='$I_5 / I_2$', color='xkcd:mustard')
    axA.plot(-b1range, jX, label='$J_x$', color='xkcd:olive')
    axA.axhline(3, color='xkcd:olive', linewidth=0.5, linestyle='dotted')
    axA.plot(-b1range, F, label='$F$', color='xkcd:royal blue')
    axB.plot(-b1range, G, label='$G$', color='black')
    axB.plot(-b1range, maxM, label='max $m$', color='xkcd:red')

    for a in (axA, axB):    
        a.set(xlim=(0,2.5), ylim=(0,10), xlabel='$b_1$')

if __name__ == '__main__':
    from sys import argv

    # if mode=='JaCoW':
    columnWidth = 3.28
    doubleWidth = 6.68
    # elif mode=='PRAB':
    #    pass


    for filename in argv[1:]:
        if filename in ("islands.pdf","Islands.pdf"):
            print("overview of 3 stability islands in k0,k1 space")
            fig, ax = subplots(figsize=(columnWidth, 1.2 * columnWidth))
            fks = FloquetKSpace(arange(-1, 4.01, 0.05), arange(0, 8.01, 0.05))
            fks.solveCxy()
            fks.plotStability(ax)
            [ax.spines[dr].set_color(None) for dr in ('top', 'right')]
            ax.set_xlim((-1, 4.01))
            fig.subplots_adjust(top=0.965, bottom=0.115, left=0.115, right=0.97)

            if filename[0]=='i':
                print('truncate to k1 <= 4')
                fig.set_size_inches(columnWidth, 0.75*columnWidth)
                ax.set_ylim((0,4))
                fig.subplots_adjust(top=0.969, bottom=0.177, left=0.127, right=0.968)
            else:
                fig.subplots_adjust(top=0.98, bottom=0.115, left=0.125, right=0.97)
            saveFig(fig, filename)

        elif filename == "necktie.pdf":
            print("zoom-in of the necktie island in k0,k1 space")
            fig, ax = subplots(figsize=(columnWidth, 0.9*columnWidth))
            fks = FloquetKSpace(arange(-0.34, 0.35, 0.02), arange(0, 0.93, 0.02))
            fks.solveCxy()
            fks.plotStability(ax, tuneLevels=arange(.1, .45, .1))
            [ax.spines[dr].set_color(None) for dr in ('top', 'right')]
            ax.yaxis.set_ticks(arange(10)/10)
            fig.subplots_adjust(top=0.98, bottom=0.14, left=0.16, right=0.975)
            saveFig(fig, filename)

        elif filename == "chroma.pdf":
            tm = generateMapIfNotExisting()
            fig, ax = subplots(figsize=(columnWidth, 0.9*columnWidth))
            grayDiagram(ax, tm, tm.chroma[:, :, 0], arange(-2.5, 0.1, 0.5), fmt='%.1f', grayDiv=5, grayMin=-1)
            ax.plot((0, 0.5), (0, 0.5), color='black', linewidth=0.3, linestyle='dashed')
            ax.set_aspect('equal')
            fig.subplots_adjust(top=0.985, bottom=0.13, left=0.157, right=0.961)
            saveFig(fig, filename)

        elif filename == "F.pdf":
            print("some figures of merit for F-optimized cells (excluding sextupole strengths)")
            tm = generateMapIfNotExisting()
            fig, ax = subplots(1, 4, figsize=(doubleWidth, 0.8*columnWidth), sharex=True, sharey=True)
            print('subplot 0: b1')
            grayDiagram(ax[0], tm, tm.mapF.b1, arange(-1.4, -1.0, 0.1), fmt='%.1f', grayDiv=5)
            print('subplot 1: F')
            plotF(ax[1], tm, tm.mapF.fval)            
            print('subplot 2: '+tm.mapF.atNames[0])
            plotJx(ax[2], tm, tm.mapF.atArray[:, :, 0])
            print('subplot 3: '+tm.mapF.atNames[2])
            grayDiagram(ax[3], tm, 1e3*tm.mapF.atArray[:, :, 2], # (-0.1, -1e-2, -1e-3, -1e-4, 1e-4, 1e-3, 1e-2, 0.1), fmt='%.1e', 
		arange(-100, 101, 25), fmt='%i',
                        faceLims=((-5, 5),), faceColors=('#ffcccc',),
                        grayDiv=5, grayMax=25, grayMin=-25)
            ax[0].set_xlim((0.2, 0.5))
            for a in ax:
                setp(a.get_xticklabels()[0], visible=False)

            fig.subplots_adjust(top=0.97, bottom=0.15, left=0.07, right=0.98, wspace=0.1)
            saveFig(fig, filename)

        elif filename in ("example.pdf", "Example.pdf"):
            fc = FourierCell()
            b1 = -1.11
            if False: # k vals from tunes
                tunes = fc.tuneTo(0.45, 0.35)
                k = fc.k
            else: # tunes from k vals
                k = array([0.04801, 0.85536])
                tunes = fc.setKxy(k)
            print("k_0 = %.5f, k_1 = %.5f, b1 = %.4f" % (*fc.k, b1))
            print("nu_x = %.8f, nu_y=%.8f" % tunes)
            fc.setB(array([b1])) #[b1,1.3]))
            print("F = %.4f, Jx = %.4f, xi_x = %.4f, xi_y = %.4f, alpha=%.6f" % 
              (fc.gr.F(),fc.gr.jX(),*fc.gr.naturalChroma(),fc.gr.momComp()))
            print("maxM = %.4f, G = %.4f" % (fc.maxM(), fc.G()))

            if filename[0]=='e':
                fig, ax = subplots(figsize=(0.5*columnWidth,0.68*columnWidth))
                sVar = 's'
                plotMultipoles(ax, fc, sVar=sVar)
                etaLabel = r'$\eta(s)$'
                betaLabel = r'$\beta_%c(s)$'
            else:
                fig, (ax, axB) = subplots(2,1,figsize=(columnWidth, columnWidth), sharex=True)
                sVar = 'u'
                plotMultipoles(axB, fc, sVar=sVar, sextupoles=True)
                etaLabel = r'$\tilde\eta(u)$'
                betaLabel = r'$\tilde\beta_%c(u)$'
                
            ax.plot(2*fc.gr.sL, fc.gr.betaX, label=betaLabel % 'x', color='xkcd:navy blue')
            ax.plot(2*fc.gr.sL, fc.gr.betaY, label=betaLabel % 'y', color='0.5')
            ax.plot(2*fc.gr.sL, fc.gr.eta, label=etaLabel, color='red')

            if filename[0]=='e':
                fig.subplots_adjust(top=0.999, bottom=0.19, left=0.164, right=0.941)
                ax.set(ylim=(-2,17))
                ax.legend(prop={'size': 8})
                [ax.spines[dr].set_color(None) for dr in ('top', 'right')]
            else:
                fig.subplots_adjust(top=0.99, bottom=0.13, left=0.1, right=0.96)
                ax.set(ylim=(0,17))
                for a in (ax, axB):
                    [a.spines[dr].set_color(None) for dr in ('top', 'right')]
                    a.legend(ncol=2)
            saveFig(fig, filename)

        elif filename == "sextuVals.pdf":
            print("F-optimized sextupole strength coefficients m0, m1")
            tm = generateMapIfNotExisting()
            fig, ax = subplots(1, 2, figsize=(columnWidth, 0.8*columnWidth), sharex=True, sharey=True)
            print('subplot 0: '+tm.mapF.atNames[3])
            grayDiagram(ax[0], tm, tm.mapF.atArray[:, :, 3], arange(-1.6, 0.1, 0.2), fmt='%.1f', grayDiv=2)
            print('subplot 1: '+tm.mapF.atNames[4])
            grayDiagram(ax[1], tm, tm.mapF.atArray[:, :, 4], arange(0, 1.4, 0.2), fmt='%.1f', grayDiv=2)
            ax[0].set_xlim((0.2, 0.5))
            for a in ax:
                setp(a.get_xticklabels()[0], visible=False)

            fig.subplots_adjust(top=0.97, bottom=0.15, left=0.155, right=0.97, wspace=0.1)
            saveFig(fig, filename)

        elif filename == "G.pdf":
            print("some figures of merit for G-optimized cells (excluding sextupole strengths)")
            tm = generateMapIfNotExisting()
            fig, ax = subplots(1, 4, figsize=(doubleWidth, 0.8*columnWidth), sharex=True, sharey=True)
            # print('subplot 0: b1')
            # grayDiagram(ax[0], tm, tm.mapG.b1, arange(-1.2, -0.1, 0.1), fmt='%.2f', grayDiv=2)
            print('subplot 0: G')
            grayDiagram(ax[0], tm, tm.mapG.fval, arange(6), fmt='%i', grayMax=2, grayDiv=10)
            print('subplot 1: '+tm.mapG.atNames[1])
            plotF(ax[1], tm, tm.mapG.atArray[:, :, 1])
            print('subplot 2: '+tm.mapG.atNames[0])
            plotJx(ax[2], tm, tm.mapG.atArray[:, :, 0])
            print('subplot 3: '+tm.mapG.atNames[2])
            grayDiagram(ax[3], tm, 1e3*tm.mapG.atArray[:, :, 2],
		arange(-100, 101, 25), fmt='%i', grayDiv=5, grayMax=25, grayMin=-25)
            # ax[0].set_xlim((0.2, 0.5))
            for a in ax:
                setp(a.get_xticklabels()[0], visible=False)

            fig.subplots_adjust(top=0.97, bottom=0.15, left=0.07, right=0.98, wspace=0.1)
            saveFig(fig, filename)

        elif filename == "b1scan.pdf":
            fig, ax = subplots(2, 2, figsize=(doubleWidth, 1.5*columnWidth), sharex=True, sharey='row')
            b1scan(ax[0,0], ax[1,0], nuX=0.15)
            b1scan(ax[0,1], ax[1,1], nuX=0.45)
            for a in ax.flat:
                [a.spines[dr].set_color(None) for dr in ('top', 'right')]
                a.label_outer()
                a.legend()
            saveFig(fig, filename, tight=True)

        elif filename == "extendedExample.pdf":
            fig, ax = subplots(figsize=(columnWidth,columnWidth))
            fc = FourierCell(mSize=3)
            fc.tuneTo(0.45,0.35)
            fc.setB()
            fc.sextuVals()
            plotMagnetStrengths

        else:
            print("unrecognized filename")
            tm = generateMapIfNotExisting()
            fig, ax = subplots()
            ax.imshow(tm.mapF.fval, origin='lower')
            saveFig(fig)

