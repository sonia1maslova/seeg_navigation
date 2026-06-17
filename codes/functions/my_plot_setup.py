import matplotlib.pyplot as plt

def set_axis (ax, axis_linewidth):
    ax.tick_params(axis='both', which='both',
                   direction='out', length=4, width=1.2,
                   bottom=True, top=False, left=True, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_linewidth(axis_linewidth)
    ax.spines['left'].set_linewidth(axis_linewidth)
    return ax


def set_plot(nrows = 1, ncols = 1, figsize=(8, 6), dpi=100,
                  fontfamily='Arial', fontsize=12, axis_linewidth=1.8,
                  label_bold=True):
    plt.rcParams.update({
        'font.family': fontfamily,
        'font.size': fontsize,
        'axes.labelsize': fontsize,
        'axes.titlesize': fontsize + 2,
        'xtick.labelsize': fontsize - 2,
        'ytick.labelsize': fontsize - 2,
        'legend.fontsize': fontsize - 1,
        'axes.labelweight': 'bold' if label_bold else 'normal',
        'axes.titleweight': 'bold' if label_bold else 'normal'
    })

    fig, ax = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi)
    if ((nrows >1) & (ncols> 1)):
        for irow in range(nrows):
            for icol in range(ncols):
                sub_ax =  ax[irow, icol]
                ax[irow, icol] = set_axis(sub_ax, axis_linewidth=axis_linewidth)
    elif ((nrows == 1) & (ncols== 1)):
        ax = set_axis(ax, axis_linewidth=axis_linewidth)
    else:
        for idx, sub_ax in enumerate(ax):
            ax[idx] = set_axis(sub_ax, axis_linewidth=axis_linewidth)

    return fig, ax