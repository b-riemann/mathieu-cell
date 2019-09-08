"""Figure generator for Mathieu cells"""
from matplotlib.pyplot import subplots

def stabilityDiagram(ax):
    """print stability diagram"""

if __name__ == '__main__':
    from sys import argv

    for mode in argv[:1]:
        fig, ax = subplots()
        if mode == 'stability':
            stabilityDiagram(ax)
            fig.savefig('stability.pdf')
        fig.close()
