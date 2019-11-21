"""Figure generator for Mathieu cells"""
from numpy import arange, squeeze
from tools import subplots, saveFig
from floquetCell import FloquetKSpace
from fourierCell import TuneMap
import os

# to be adapted to fourierCell.TuneMap (before: floquetCell.TuneMap)
def generateMapIfNotExisting(npz_filename="tmap.npz", npz_cubefile="cubes.npz"):
    tm = TuneMap() # b1Range=-arange(0, 1.5, 0.02)) # -arange(0, 1.5, 0.05)  # arange(0, 2, 0.1), arange(-1.4, 0.01, 0.05)
    # if os.path.isfile(npz_filename):
    #     tm.loadKs(npz_filename)
    # else:
    #     print("computing k values for tune map. this only needs to be done once.")
    #     tm.generateKs()
    #     tm.saveKs(npz_filename)

    # if os.path.isfile(npz_cubefile):
    #     tm.loadCubes(npz_cubefile)
    # else:
    #     print("computing objective functions on cube. this only needs to be done once.")
    #     tm.generateCubes()
    #     tm.saveCubes(npz_cubefile)
    # tm.processCubes()
    tm.make()
    return tm

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
        
        elif filename == "tunemap-check.pdf":
            print("implement this using the new TuneMap class in fourierCell.pyx")
            tm = generateMapIfNotExisting()
            fig, ax = subplots(figsize=(4, 4))
            tm.grayDiagram(ax, squeeze(tm.cubes['tuneX']), arange(0.025, 0.49, 0.025), grayDiv=None)
            tm.grayDiagram(ax, squeeze(tm.cubes['tuneY']), arange(0.025, 0.49, 0.025), grayDiv=None)
            ax.set_aspect('equal')
            saveFig(fig, filename, tight=True)

        else:
            print("unrecognized filename")
