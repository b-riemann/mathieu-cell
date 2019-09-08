"""Figure generator for Mathieu cells"""
from numpy import arange
from matplotlib.pyplot import subplots
from floquetCell import FloquetKSpace

def stabilityDiagram(ax):
    """print stability diagram"""
    fks = FloquetKSpace(arange(-1, 4.01, 0.05), arange(0, 8.01, 0.05))
    fks.solveCxy()
    fks.plotStability(ax)

if __name__ == '__main__':
    from sys import argv

    for mode in argv[:1]:
        if mode == 'stability':
            fig, ax = subplots()
            stabilityDiagram(ax)
            fig.savefig('stability.pdf')
            fig.close()
