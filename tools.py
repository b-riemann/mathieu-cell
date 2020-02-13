from matplotlib.pyplot import subplots, show, close, clabel
from matplotlib import rcParams
from numpy import arange, savez, load, mod

rcParams['contour.negative_linestyle'] = 'solid'
rcParams.update({'font.family': 'serif', 'font.size': 9})

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


def plot_multp(ax, u, b, k, m=None, sVar='s', dipoleColor='black', quadColor='xkcd:ocean blue', sextColor='xkcd:mustard'):
    ax.plot(u, b, label='$b(%c)$' % sVar, linewidth=0.7, color=dipoleColor)
    ax.plot(u, k, label='$k(%c)$' % sVar, linewidth=0.7, color=quadColor)
    if m is not None:
        ax.plot(u, m, label='$m(%c)$' % sVar, linewidth=0.7, color=sextColor)
    [ax.spines[dr].set_color(None) for dr in ('top', 'right')]


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


