"""Figure generator for Mathieu cells"""
from matplotlib import rc
from matplotlib.pyplot import setp
from numpy import arange, squeeze, array, ceil, deg2rad, rad2deg, diff, empty_like, amax, absolute, NaN, pi, mean, ones_like, outer, sin, sqrt
from scipy.optimize import minimize as minix, minimize_scalar as minix_scalar
from tools import subplots, saveFig, centaur
from floquetCell import FloquetKSpace
from fourierCell import TuneMap, FourierCell, Grapher
from time import time
import os

# rc('text', usetex=True)
rc('font', family='serif', size=8)
# rc('contour', negative_linestyle='solid')

def generateMapIfNotExisting():
    tm = TuneMap()
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
                faceLims=((-10, 0), (3, 10000)), faceColors=('#cccccc', '#cccccc'))


def plot_multp(ax, u, b, k, m=None, sVar='s', dipoleColor='black', quadColor='xkcd:ocean blue', sextColor='xkcd:mustard'):
    ax.plot(u, b, label='$b(%c)$' % sVar, linewidth=0.7, color=dipoleColor)
    ax.plot(u, k, label='$k(%c)$' % sVar, linewidth=0.7, color=quadColor)
    if m is not None:
        ax.plot(u, m, label='$m(%c)$' % sVar, linewidth=0.7, color=sextColor)
    [ax.spines[dr].set_color(None) for dr in ('top', 'right')]


def plotMultipoles(ax, gr : Grapher, sVar='s', sextupoles=False):
    plot_multp(ax, 2*gr.sL, gr.b, gr.k, gr.mSext if sextupoles else None, sVar=sVar)
    ax.axhline(0, color='black', linestyle='dashed', linewidth=0.5)
    ax.set(xlim=(0,1), xlabel=r'%c / $2\pi$' % sVar)


def obig(b2, b1, fc : FourierCell):
    fc.setB(array([b1,b2]))
    return fc.G()


def b1scan(axA, axB, nuX=0.45, nuY=0.35, minim=False, b1range=-arange(0,2.5,0.02)):
    fc = FourierCell()
    fcTri = FourierCell(mSize=3)
    tunes = fc.tuneTo(nuX, nuY)
    tunes = fcTri.tuneTo(nuX, nuY)
    print("k_0 = %.5f, k_1 = %.5f" % tuple(fc.k))
    F = empty_like(b1range)
    jX = empty_like(b1range)
    i5i2 = empty_like(b1range)
    i1 = empty_like(b1range)
    G = empty_like(b1range)
    maxM = empty_like(b1range)
    Gtri = empty_like(b1range)
    maxMtri = empty_like(b1range)
    for n, b1 in enumerate(b1range):
        if minim:
            # result = minix(obig, 0, args=(b1,fc))
            result = minix_scalar(obig, bounds=(-2,2), args=(b1,fc))
            fc.setB(array([b1,result.x]))
            fcTri.setB(array([b1,result.x]))
        else:
            fc.setB(array([b1]))
            fcTri.setB(array([b1]))
        F[n] = fc.gr.F()
        G[n] = fc.G()
        jX[n] = fc.gr.jX()
        i5i2[n] = fc.gr.i5sum() / fc.gr.i2sum()
        i1[n] = fc.gr.i1()
        maxM[n] = amax(absolute(fc.gr.mSext))
        Gtri[n] = fcTri.G()
        maxMtri[n] = amax(absolute(fcTri.gr.mSext))
        # print('b1 = %.2f : chrom.nat.x=%.4f,y=%.4f, chrom.full.x=%.4f,y=%.4f' % (b1, *fc.gr.naturalChroma(), *fc.gr.fullChroma()))
    F[jX <= 0] = NaN
    G[jX <= 0] = NaN
    Gtri[jX <= 0] = NaN

    # axA.plot(-b1range, i5i2, label='$I_5 / I_2$', linewidth=0.7, color='xkcd:mustard')
    axA.plot(-b1range, i1, label='$I_1$', linewidth=0.7, color='xkcd:mustard')
    axA.plot(-b1range, jX, label='$J_x$', linewidth=0.7, color='xkcd:olive')
    axA.axhline(3, color='black', linewidth=0.5, linestyle='dotted')
    axA.axhline(0, color='black', linewidth=0.5, linestyle='dotted')
    axA.plot(-b1range, F, label='$F$', linewidth=0.7, color='xkcd:royal blue')
    axB.plot(-b1range, G, label='$G$', linewidth=0.7, color='xkcd:mustard')
    axB.plot(-b1range, Gtri, label='$G_{tri}$', linewidth=0.7, color='magenta')
    axB.plot(-b1range, maxM, label='max $m$', linewidth=0.7, color='xkcd:red')
    axB.plot(-b1range, maxMtri, label='max $m_{tri}$', linewidth=0.7, color='0.3')

    axA.set(xlim=(0,2.5), ylim=(-1,9.5), xlabel='$b_1$')
    axB.set(xlim=(0,2.5), ylim=(0,6), xlabel='$b_1$')
    for a in (axA, axB):    
        [a.spines[dr].set_color(None) for dr in ('top', 'right')]


def multicolor_ylabel(ax, list_of_strings, list_of_colors, anchorpad=0, **kw):
    """from https://stackoverflow.com/questions/33159134/matplotlib-y-axis-label-with-multiple-colors
    this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    boxes = [TextArea(text, textprops=dict(color=color, ha='left', va='bottom', rotation=90, **kw))
             for text, color in zip(list_of_strings[::-1], list_of_colors[::-1])]
    ybox = VPacker(children=boxes, align="center", pad=0, sep=10)
    anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.2, 0.3),
                                      bbox_transform=ax.transAxes, borderpad=0.)
    ax.add_artist(anchored_ybox)


def opticsPlot(axBeta, s, betaX, betaY, etaX, betaXcolor='xkcd:royal blue', betaYcolor='0.4', etaColor='red', prefix='', etaLim=2.5):
    axEta = axBeta.twinx() 
    axBeta.plot(s, betaX, color=betaXcolor, label='x', linewidth=1)
    axBeta.plot(s, betaY, color=betaYcolor, label='y', linewidth=1)
    multicolor_ylabel(axBeta, (r'$%s\beta_x$' % prefix, r'$%s\beta_y$' % prefix), (betaXcolor, betaYcolor))
    axEta.plot(s, etaX, color=etaColor, linewidth=1)
    axEta.set(ylim=(0, etaLim), yticks=arange(0, etaLim + 1))
    axEta.set_ylabel(r'$%s\eta_x$' % prefix, color=etaColor)

    axBeta.spines['top'].set_color(None)
    axEta.spines['top'].set_color(None)
    axBeta.spines['right'].set_color(None)
    axEta.spines['left'].set_color(None)


def opaExport(filename, gr : Grapher, cellLength, energyGeV=2.4, curvature=deg2rad(5)/2.4):  # s, b, kArr, mArr):
    s = gr.sL * cellLength
    scaler = cellLength / pi

    elemLength = mean(diff(s))
    bendAngles = rad2deg(gr.b * elemLength) * curvature
    kArr = gr.k / scaler**2
    intSextuStrength = gr.mSext * (elemLength / 2) / curvature / scaler**4
    lineElems = list()

    with open(filename, 'w', encoding='cp1252') as f:
        f.write('energy = %.6f;\n{--- table of elements ---}\n' % energyGeV)
        for n, (t, k, mL) in enumerate(zip(bendAngles, kArr, intSextuStrength)):  # fixit later
            f.write('thick%i : combined, l=%.6f, t=%.6f, k=%.6f, ax=10.00, ay=10.00;\n' % (n, elemLength, t, k))
            f.write('sextu%i : sextupole, k=%.6f, ax=10.00, ay=10.00;\n' % (n, mL))
            lineElems.extend(['thick%i' % n, 'sextu%i' % n])

        f.write('cell : %s;\n' % ','.join(lineElems))
    print('sliced lattice exported to %s' % filename)


def poleSheets(gr : Grapher, weighted=False, phiSteps=16):
    phi= arange(phiSteps)/phiSteps * 2*pi
    weights = (1,2,3) if weighted else (1,1,1)
    return weights[0]*outer(sin(phi), gr.b), weights[1]*outer(sin(2*phi), gr.k), weights[2]*outer(sin(3*phi), gr.mSext)


def showSLSparameters(maxM=None):
    print("SLS parameters:")
    cellAngle = deg2rad(5.0)
    cellLength = 2.165 # [m]
    iinvRho = cellLength / cellAngle # [m]
    bRho = 8.0
    characteristicB = bRho/iinvRho
    radius = 0.01
    characteristicLength = pi*sqrt(radius*iinvRho)
    print(r"b \rho = %.1f T m" % bRho)
    print(r"1 / < 1 / \rho > = %.2f m" % iinvRho)
    print(r"B_c = %.4f T" % characteristicB)
    print(r"L_c = %.4f m (R=%.3f m) " % (characteristicLength, radius))

    maxMu = 650.0
    peakFld = maxMu * bRho * radius**2
    print("max |mu| = %.1f, peak at %.1f mm = %.3f T" % (maxMu, radius*1000, peakFld))

    if maxM is not None:
        optLength = pi * (maxM * iinvRho / maxMu)**0.25
        print("max |m| = %.4f: optimal length = %.3f m" % (maxM, optLength))
    return iinvRho, characteristicB, characteristicLength


def poleTipVals(ax, gr : Grapher, LcL_range=arange(0,1.81,0.01), weighted=False):
    sheets = poleSheets(gr, weighted=False)
    sheets_w = poleSheets(gr, weighted=True)

    max_m = amax(absolute(gr.mSext))
    print('max m (by poletip) = %.4f' % max_m)
    BrBc_m = empty_like(LcL_range) # only sextupole
    BrBc_bkm = empty_like(LcL_range)
    BrBc_bkm_w = empty_like(LcL_range)
    for n, LcL in enumerate(LcL_range):  # i know this can be done without loop, but its still fast enough
        BrBc_m[n] = max_m*LcL**4
        BrBc_surf = sheets[0] + sheets[1]*LcL**2 + sheets[2]*LcL**4  
        BrBc_bkm[n] = amax(absolute(BrBc_surf)) 
        BrBc_surf_w = sheets_w[0] + sheets_w[1]*LcL**2 + sheets_w[2]*LcL**4  
        BrBc_bkm_w[n] = amax(absolute(BrBc_surf_w)) 
    
    _, characteristicB, characteristicLength = showSLSparameters()
    print(r"B_r / B_c = %.3f, max B_dipole(b1) = %.2f T" % (BrBc_bkm[0], BrBc_bkm[0]*characteristicB))

    maxField = 2.0  # T
    bMax = maxField / characteristicB
    print(r"max B = %.1f T,  max b = %.4f, max |b1| = %.4f" % (maxField, bMax, (bMax-1)/2))
    ax.axhline(bMax, linewidth=0.5, color='black', linestyle='dotted')

    LcL_ref = 0.855

    ax.plot(LcL_range, BrBc_bkm_w, color='red')
    ax.plot(LcL_range, BrBc_bkm, color='0.6')
    ax.plot(LcL_range, BrBc_m, linewidth=0.7, color='xkcd:mustard')
    ax.axvline(LcL_ref, linewidth=0.5, color='black', linestyle='dotted')
    ax.set(xlabel=r'$L_c$  / $L$', ylabel=r'max $|B_r|$  / $B_c$', 
           xlim=(0,LcL_range[-1]), ylim=(0,25), yticks=arange(0,26,5))
    [ax.spines[dr].set_color(None) for dr in ('top', 'right')]

    lOPA = ceil(1000*characteristicLength/LcL_ref )/1000
    print('L_opa = %.3f m' % lOPA)
    return LcL_ref, lOPA


def poleTipContribs(ax, gr : Grapher, LcL, characteristicB, lOPA):
    s = lOPA * gr.sL

    sheets = poleSheets(gr, weighted=False)
    BrBc_surf = sheets[0] + sheets[1]*LcL**2 + sheets[2]*LcL**4  
    sheets = poleSheets(gr, weighted=True)
    BrBc_surf_weighted = sheets[0] + sheets[1]*LcL**2 + sheets[2]*LcL**4  


    # we want to show actual pole-tip field of multipoles w.o. weighting:
    plot_multp(ax, s, absolute(gr.b)*characteristicB,
                      absolute(gr.k)*characteristicB*LcL**2, 
                      absolute(gr.mSext)*characteristicB*LcL**4)
    ax.plot(s, characteristicB*amax(absolute(BrBc_surf), axis=0), color='0.6')
    ax.plot(s, characteristicB*amax(absolute(BrBc_surf_weighted), axis=0), color='red')
    ax.axhline(2.0, linewidth=0.5, color='black', linestyle='dotted')
    ax.set(xlim=(0,lOPA), xlabel=r'$s$ [m]', ylim=(0,2.02), ylabel=r'max $|B|$ [T]')


if __name__ == '__main__':
    from sys import argv

    columnWidth = 3.28
    doubleWidth = 7.0  # 6.68

    exampleA = (0.45,0.35)
    exampleB = (0.15,0.35)

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
            grayDiagram(ax[0], tm, -tm.mapF.b1, arange(1.0, 1.4, 0.1), fmt='%.1f', grayDiv=5)
            print('subplot 1: F')
            plotF(ax[1], tm, tm.mapF.fval)            
            print('subplot 2: '+tm.mapF.atNames[0])
            plotJx(ax[2], tm, tm.mapF.atArray[:, :, 0])
            print('subplot 3: '+tm.mapF.atNames[2])
            grayDiagram(ax[3], tm, tm.mapF.atArray[:, :, 2], # (-0.1, -1e-2, -1e-3, -1e-4, 1e-4, 1e-3, 1e-2, 0.1), fmt='%.1e', 
		arange(-1, 1.1, 0.25), fmt='%.2f',
                        faceLims=((-0.05, 0.05),), faceColors=('#ffcccc',),
                        grayDiv=5, grayMax=0.25, grayMin=-0.25)
            ax[0].set_xlim((0.2, 0.5))
            for a in ax:
                setp(a.get_xticklabels()[0], visible=False)

            fig.subplots_adjust(top=0.97, bottom=0.15, left=0.075, right=0.985, wspace=0.045)
            saveFig(fig, filename)

        elif filename in ("example.pdf", "Example.pdf"):
            fc = FourierCell()
            b1 = -1.11
            if False: # k vals from tunes
                tunes = fc.tuneTo(*exampleA)
                k = fc.k
            else: # tunes from k vals
                k = array([0.04801, 0.85536])
                tunes = fc.setKxy(k)
            print("k_0 = %.5f, k_1 = %.5f, b1 = %.4f" % (*fc.k, -b1))
            print("nu_x = %.8f, nu_y=%.8f" % tunes)
            fc.setB(array([b1])) #[b1,1.3]))
            print("F = %.4f, Jx = %.4f, xi_x = %.4f, xi_y = %.4f, I1=%.6f" % 
              (fc.gr.F(),fc.gr.jX(),*fc.gr.naturalChroma(),fc.gr.i1()))
            print("m0=%.4f, m1=%.4f, G = %.4f" % (*fc.sextuVals(), fc.G()))

            showSLSparameters( fc.maxM() )

            if filename[0]=='e':
                fig, ax = subplots(figsize=(0.5*columnWidth,0.68*columnWidth))
                sVar = 's'
                plotMultipoles(ax, fc.gr, sVar=sVar)
                betaLabel = r'$\beta_%c(s)$'
                ax.plot(2*fc.gr.sL, fc.gr.betaX, label=betaLabel % 'x', color='xkcd:royal blue')
                ax.plot(2*fc.gr.sL, fc.gr.betaY, label=betaLabel % 'y', color='0.4')
                ax.plot(2*fc.gr.sL, fc.gr.eta, label=r'$\eta(s)$', color='red')

            else:
                fig, (ax, axB) = subplots(2,1,figsize=(columnWidth, columnWidth), sharex=True)
                sVar = 'u'
                plotMultipoles(axB, fc.gr, sVar=sVar, sextupoles=True)
                opticsPlot(ax, 2*fc.gr.sL, fc.gr.betaX, fc.gr.betaY, fc.gr.eta, prefix=r'\tilde') 

            if filename[0]=='e':
                fig.subplots_adjust(top=0.999, bottom=0.19, left=0.164, right=0.941)
                ax.set(ylim=(-2,17))
                ax.legend(prop={'size': 8})
                [ax.spines[dr].set_color(None) for dr in ('top', 'right')]
            else:
                fig.subplots_adjust(top=0.99, bottom=0.13, left=0.15, right=0.87)
                ax.set(ylim=(0,17))

            print("nat. chroma %.4f %.4f, full chroma %.1e %.1e" % (*fc.gr.naturalChroma(), *fc.gr.fullChroma()))
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
            # grayDiagram(ax[0], tm, -tm.mapG.b1, arange(0.1, 1.2, 0.1), fmt='%.2f', grayDiv=2)
            print('subplot 0: G')
            grayDiagram(ax[0], tm, tm.mapG.fval, arange(6), fmt='%i', grayMax=2, grayDiv=10)
            print('subplot 1: '+tm.mapG.atNames[1])
            plotF(ax[1], tm, tm.mapG.atArray[:, :, 1])
            print('subplot 2: '+tm.mapG.atNames[0])
            plotJx(ax[2], tm, tm.mapG.atArray[:, :, 0])
            print('subplot 3: '+tm.mapG.atNames[2])
            grayDiagram(ax[3], tm, tm.mapG.atArray[:, :, 2],
		arange(-1, 1.1, 0.25), fmt='%.2f', grayDiv=5, grayMax=0.25, grayMin=-0.25)
            # ax[0].set_xlim((0.2, 0.5))
            for a in ax:
                setp(a.get_xticklabels()[0], visible=False)

            fig.subplots_adjust(top=0.97, bottom=0.15, left=0.075, right=0.985, wspace=0.045)
            saveFig(fig, filename)

        elif filename.startswith("b1scan"):
            fig, ag = subplots(figsize=(columnWidth,0.7*columnWidth))
            fih, ah = subplots(figsize=(columnWidth,0.7*columnWidth))
            if filename.startswith("b1scanA"):
                chara = 'A'
                b1scan(ag, ah, *exampleA) #, b1range=-arange(0.9,1.3,0.01))
                [a.axvline(1.11, color='black', linewidth=0.5, linestyle='dotted') for a in (ag, ah)]
                # ag.set(xlim=(0.9,1.3),ylim=(-1,15))
            else:
                chara = 'B'
                b1scan(ag, ah, *exampleB)
            fig.subplots_adjust(top=0.96, bottom=0.19, left=0.07, right=0.96)
            saveFig(fig, "b1scan%c-int.pdf" % chara)
            fih.subplots_adjust(top=0.96, bottom=0.19, left=0.07, right=0.96)
            saveFig(fih, "b1scan%c-mG.pdf" % chara)
            
        elif filename in ("extendedExample.pdf", "poleTip.pdf", "example.opa"):
            fc = FourierCell(mSize=3)
            b1 = -1.11
            
            tunes = fc.tuneTo(0.45, 0.35)
            fc.setB(array([b1]))
            fc.sextuVals()
            _, characteristicB, _ = showSLSparameters( fc.maxM() )
            
            fig, ax = subplots(figsize=(columnWidth,0.6*columnWidth))
            plotMultipoles(ax, fc.gr, 'u', sextupoles=True)
            print("nat. chroma %.4f %.4f, full chroma %.1e %.1e" % (*fc.gr.naturalChroma(), *fc.gr.fullChroma()))
            fig.subplots_adjust(top=0.98, bottom=0.225, left=0.15, right=0.87)
            saveFig(fig, "extendedExample.pdf")

            fig, ax = subplots(figsize=(columnWidth,0.7*columnWidth))
            LcL_ref, lOPA = poleTipVals(ax, fc.gr, weighted=True)
            fig.subplots_adjust(top=.98,bottom=.16,left=.16,right=.98)
            saveFig(fig, "poleTip.pdf") 

            opaExport("example.opa", fc.gr, lOPA)

            fig, ax = subplots(figsize=(columnWidth,0.9*columnWidth))
            poleTipContribs(ax, fc.gr, LcL_ref, characteristicB, lOPA)
            fig.subplots_adjust(top=.98,bottom=.16,left=.18,right=.98)
            saveFig(fig, "poleTipContribs.pdf")


        elif filename == "H.pdf":
            print("some figures of merit for H-optimized cells (excluding sextupole strengths)")
            tm = generateMapIfNotExisting()
            fig, ax = subplots(1, 4, figsize=(doubleWidth, 0.8*columnWidth), sharex=True, sharey=True)
            # print('subplot 0: b1')
            # grayDiagram(ax[0], tm, tm.mapH.b1, arange(-1.2, -0.1, 0.1), fmt='%.2f', grayDiv=2)
            print('subplot 0: G')
            grayDiagram(ax[0], tm, tm.mapH.fval, arange(6), fmt='%i', grayMax=2, grayDiv=10)
            print('subplot 1: '+tm.mapH.atNames[1])
            plotF(ax[1], tm, tm.mapH.atArray[:, :, 1])
            print('subplot 2: '+tm.mapH.atNames[0])
            plotJx(ax[2], tm, tm.mapH.atArray[:, :, 0])
            print('subplot 3: '+tm.mapH.atNames[2])
            grayDiagram(ax[3], tm, tm.mapH.atArray[:, :, 2],
		arange(-1, 1.1, 0.25), fmt='%.2f', grayDiv=5, grayMax=0.25, grayMin=-0.25)
            # ax[0].set_xlim((0.2, 0.5))
            for a in ax:
                setp(a.get_xticklabels()[0], visible=False)

            fig.subplots_adjust(top=0.97, bottom=0.15, left=0.075, right=0.985, wspace=0.045)
            saveFig(fig, filename)

        else:
            print("unrecognized filename")
            tm = generateMapIfNotExisting()
            fig, ax = subplots()
            ax.imshow(tm.mapF.fval, origin='lower')
            saveFig(fig)

