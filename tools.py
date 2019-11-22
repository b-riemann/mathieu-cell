from matplotlib.pyplot import subplots, show, close, clabel
from matplotlib import rcParams
from numpy import arange, savez, load, mod

rcParams['contour.negative_linestyle'] = 'solid'

def saveFig(fig, filename=None, tight=False, **kwargs):
    if tight:
        fig.tight_layout()
    if filename is not None:
        fullname = filename
        fig.savefig(fullname, **kwargs)
        print('figure saved to %s' % fullname)
    else:
        show()
    close(fig)


def saveData(filename, **kwargs):
    fullname = filename  #'../data/%s' % filename
    print('data saved to %s' % fullname)
    savez(fullname, **kwargs)


def loadData(filename):
    fullname = filename  #'../data/%s' % filename
    print('loading data from %s' % fullname)
    return load(fullname)

def absq(x):
    return x.real ** 2 + x.imag ** 2

## nice contour plotting, linewidths=1?
def centaur(ax, x, y, z, levels, linewidths=0.7, colors='black',
            fmt=None, fontsize=8,
            grayDiv=None, grayMin=None, grayMax=None, grayColors='0.4', grayLinewidths=0.6,
            faceLims=None, faceColors=None, **kwargs):
    patch = ax.contour(x,y,z,levels, linewidths=linewidths, colors=colors, **kwargs)
    if fmt is not None:
        clabel(patch, fmt=fmt, inline=True, fontsize=fontsize)
    if grayDiv is not None:
        if grayMin is None:
            grayMin = levels[0]
        if grayMax is None:
            grayMax = levels[-1]
        grayLvl = arange(grayMin, grayMax, (levels[1] - levels[0]) / grayDiv)
        grayLevels = [grayLvl[n] for n in range(len(grayLvl)) if mod(n, grayDiv) != 0]
        ax.contour(x,y,z, grayLevels, linewidths=grayLinewidths, colors=grayColors)
    if faceLims is not None:
        for faceLim, faceColor in zip(faceLims, faceColors):
            ax.contourf(x,y,z, faceLim, colors=faceColor)
