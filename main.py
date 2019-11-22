"""Figure generator for Mathieu cells"""
from numpy import arange, squeeze
from tools import subplots, saveFig, centaur
from floquetCell import FloquetKSpace
from fourierCell import TuneMap
import os


# to be adapted to fourierCell.TuneMap (before: floquetCell.TuneMap)
def generateMapIfNotExisting():
    tm = TuneMap() # b1Range=-arange(0, 1.5, 0.02)) # -arange(0, 1.5, 0.05)  # arange(0, 2, 0.1), arange(-1.4, 0.01, 0.05)
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


def grayDiagram(ax, tm : TuneMap, quantity, levels=arange(11), fmt='%i', grayDiv=5, **kwargs):
    centaur(ax, tm.tuneX.flat, tm.tuneY.flat, quantity, levels, fmt=fmt, grayDiv=grayDiv, **kwargs)
    ax.set(xlabel=r'$\nu_x$', ylabel=r'$\nu_y$', xlim=(0, 0.5), ylim=(0, 0.5))
    ax.label_outer()


if __name__ == '__main__':
    from sys import argv
    
    for filename in argv[1:]:
        if filename == "islands.pdf":
            print("overview of 3 stability islands in k0,k1 space")
            fig, ax = subplots()
            fks = FloquetKSpace(arange(-1, 4.01, 0.05), arange(0, 8.01, 0.05))
            fks.solveCxy()
            fks.plotStability(ax)
            saveFig(fig, filename, tight=True)

        elif filename == "necktie.pdf":
            print("zoom-in of the necktie island in k0,k1 space")
            fig, ax = subplots(figsize=(4, 3.5))
            fks = FloquetKSpace(arange(-0.26, 0.27, 0.02), arange(0, 0.93, 0.02))
            fks.solveCxy()
            fks.plotStability(ax, tuneLevels=arange(.1, .45, .1))
            saveFig(fig, filename, tight=True)
        
        elif filename == "chroma.pdf":
            tm = generateMapIfNotExisting()
            fig, ax = subplots(figsize=(4, 4))
            grayDiagram(ax, tm, tm.chroma[:,:,0], arange(-2.5, 0.1, 0.5), fmt='%.1f', grayDiv=5, grayMin=-1)
            ax.plot((0, 0.5), (0, 0.5), color='xkcd:royal blue', linewidth=0.2)
            ax.set_aspect('equal')
            fig.subplots_adjust(top=0.98, bottom=0.12, left=0.14, right=0.96, hspace=0.2, wspace=0.2)
            saveFig(fig, filename)

        elif filename == "F.pdf":
            tm = generateMapIfNotExisting()
            fig, ax = subplots(1, 2, figsize=(6, 4), sharex=True, sharey=True)
            grayDiagram(ax[0], tm, tm.mapF.b1, arange(-1.4, -1.0, 0.1), fmt='%.1f', grayDiv=5)
            grayDiagram(ax[1], tm, tm.mapF.fval, arange(6), fmt='%i', grayDiv=5,
                        faceLims=((0.0,1.0),), faceColors=('#cccccc',))
            ax[0].set_xlim((0.25,0.5))
            fig.subplots_adjust( top=0.948, bottom=0.11, left=0.1, right=0.96, hspace=0.2, wspace=0.16)
            saveFig(fig, filename)

        else:
            print("unrecognized filename")
