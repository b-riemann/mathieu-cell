"""Figure generator for Mathieu cells"""
from matplotlib.pyplot import setp
from numpy import arange, squeeze
from tools import subplots, saveFig, centaur
from floquetCell import FloquetKSpace
from fourierCell import TuneMap
import os


# to be adapted to fourierCell.TuneMap (before: floquetCell.TuneMap)
def generateMapIfNotExisting():
    tm = TuneMap()  # b1Range=-arange(0, 1.5, 0.02)) # -arange(0, 1.5, 0.05)  # arange(0, 2, 0.1), arange(-1.4, 0.01, 0.05)
    if os.path.isfile('F.pkl'):
        tm.load()
    else:
        print("computing k values for tune map. this only needs to be done once.")
        tm.make()
    #     tm.generateKs()
    #     tm.saveKs(npz_filename)

    # if os.path.isfile(npz_cubefile):
    #     tm.loadCubes(npz_cubefile)
    # else:
    #     print("computing objective functions on cube. this only needs to be done once.")
    #     tm.generateCubes()
    #     tm.saveCubes(npz_cubefile)
    # tm.processCubes()
    return tm


def grayDiagram(ax, tm: TuneMap, quantity, levels=arange(11), fmt='%i', grayDiv=5, borderColors=('blue', 'green'),
                **kwargs):
    centaur(ax, tm.tuneX.flat, tm.tuneY.flat, quantity, levels, fmt=fmt, grayDiv=grayDiv, **kwargs)
    ax.set(xlabel=r'$\nu_x$', ylabel=r'$\nu_y$', xlim=(0, 0.5), ylim=(0, 0.5))
    ax.spines['right'].set_color(borderColors[0])
    ax.spines['top'].set_color(borderColors[1])
    ax.label_outer()


if __name__ == '__main__':
    from sys import argv

    columnWidth = 3.28  # adjusted for JaCoW
    doubleWidth = 6.68

    for filename in argv[1:]:
        if filename == "islands.pdf":
            print("overview of 3 stability islands in k0,k1 space")
            fig, ax = subplots(figsize=(columnWidth, 1.2 * columnWidth))
            fks = FloquetKSpace(arange(-1, 4.01, 0.05), arange(0, 8.01, 0.05))
            fks.solveCxy()
            fks.plotStability(ax)
            [ax.spines[dr].set_color(None) for dr in ('top', 'right')]
            ax.set_xlim((-1, 4.01))
            fig.subplots_adjust(top=0.965, bottom=0.115, left=0.115, right=0.97, hspace=0.2, wspace=0.2)
            saveFig(fig, filename)

        elif filename == "necktie.pdf":
            print("zoom-in of the necktie island in k0,k1 space")
            fig, ax = subplots(figsize=(columnWidth, columnWidth))
            fks = FloquetKSpace(arange(-0.26, 0.27, 0.02), arange(0, 0.93, 0.02))
            fks.solveCxy()
            fks.plotStability(ax, tuneLevels=arange(.1, .45, .1))
            [ax.spines[dr].set_color(None) for dr in ('top', 'right')]
            fig.subplots_adjust(top=0.97, bottom=0.13, left=0.16, right=0.975, hspace=0.2, wspace=0.2)
            saveFig(fig, filename)

        elif filename == "chroma.pdf":
            tm = generateMapIfNotExisting()
            fig, ax = subplots(figsize=(columnWidth, columnWidth))
            grayDiagram(ax, tm, tm.chroma[:, :, 0], arange(-2.5, 0.1, 0.5), fmt='%.1f', grayDiv=5, grayMin=-1)
            ax.plot((0, 0.5), (0, 0.5), color='black', linewidth=0.3, linestyle='dashed')
            ax.set_aspect('equal')
            fig.subplots_adjust(top=0.975, bottom=0.117, left=0.157, right=0.961, hspace=0.2, wspace=0.2)
            saveFig(fig, filename)

        elif filename == "F.pdf":
            tm = generateMapIfNotExisting()
            fig, ax = subplots(1, 4, figsize=(doubleWidth, 0.8*columnWidth), sharex=True, sharey=True)
            print('subplot 0: b1')
            grayDiagram(ax[0], tm, tm.mapF.b1, arange(-1.4, -1.0, 0.1), fmt='%.1f', grayDiv=5)
            print('subplot 1: F')
            grayDiagram(ax[1], tm, tm.mapF.fval, arange(6), fmt='%i', grayMax=2,
                        faceLims=((0.0, 1.0),), faceColors=('#ffffcc',))
            print('subplot 2: '+tm.mapF.atNames[0])
            grayDiagram(ax[2], tm, tm.mapF.atArray[:, :, 0], arange(0, 4, 0.5), fmt='%.1f', grayMax=2.5,
                        faceLims=((-10, 0), (3, 10)), faceColors=('#cccccc', '#cccccc'))
            print('subplot 3: '+tm.mapF.atNames[2])
            grayDiagram(ax[3], tm, tm.mapF.atArray[:, :, 2], arange(-0.5, 0.55, 0.1), fmt='%.1f',
                        faceLims=((-0.02, 0.02),), faceColors=('#ffcccc',),
                        grayDiv=5, grayMax=0.2, grayMin=-0.2)
            ax[0].set_xlim((0.2, 0.5))
            for a in ax:
                setp(a.get_xticklabels()[0], visible=False)

            fig.subplots_adjust(top=0.97, bottom=0.15, left=0.07, right=0.98, wspace=0.1)
            saveFig(fig, filename)

        else:
            print("unrecognized filename")
