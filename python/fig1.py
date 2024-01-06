"""Plot the Fig. 1 in the paper.
"""
import csv
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import PercentFormatter
import matplotlib.gridspec as gridspec
import seaborn as sns


PIC_NAME_UNIT = 'fig1.png'

FIG_SIZE_SINGLE_COLUMN = (8, 4)
LABEL_SIZE_SINGLE_COLUMN = 12
# TITLE_SIZE_SINGLE_COLUMN = 16
TITLE_SIZE_SINGLE_COLUMN = 14

FIG_SIZE_DOUBLE_COLUMN = (16, 4)
LABEL_SIZE_DOUBLE_COLUMN = 10
TITLE_SIZE_DOUBLE_COLUMN = 12

LABEL_SIZE_3D = 10
TITLE_SIZE_3D = 12

LINE_STYLE = '-'
LINE_WIDTH = 1.


COLOR_PALETTE = ["#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b",
                 "#bdcf32", "#87bc45", "#00bfa0", "#27aeef", "#6481D9"]
MARKER_STYLE = "d"
MARKER_SIZE = 4
MARKER_SIZE_BIG = 30

Y_SUBTITLE = -0.4
IDS = np.arange(1, N+1, 1)

my_cmap = sns.diverging_palette(120., 255., s=100, as_cmap=True)

def plot_scatter():
    """The scatter plot for the RIM results and the closed-form solutions.
    """

    inter_round = []
    allocation_gd = []
    allocation_fit = []
    
    with open(FILE_NAME, 'r', encoding='utf-8') as file:
        reader = csv.reader(file)

        count = 0
        for row in reader:

            count = count + 1

    x_grid = np.linspace(0, 700, 1000)
    Y_grid = np.array(x_grid)

    with np.load(FILE_NAME_PVCG_RESULT) as data:
        payments_rim = data['array1']
        x_vec_rim = data['array2']
        z_matrix_rim = data['array3']

    z_diff = (z_matrix_rim - z_matrix_closed)/ z_matrix_closed
    for i in range(N):
        z_diff[i,i] = 0.

    fig = plt.figure(figsize=(9,4))
    gs = gridspec.GridSpec(
        1, 3, 
        width_ratios=(1.,1.,1.3), wspace=0.3) 

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.set_xlabel('Closed-Form Solution', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    ax2.set_xlabel('Closed-Form Solution', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    
    ax1.set_ylabel('RIM Result', size = LABEL_SIZE_SINGLE_COLUMN)    
    ax2.set_ylabel('RIM Result', size = LABEL_SIZE_SINGLE_COLUMN)    
    
    ax1.set_aspect('equal', 'box') 
    ax2.set_aspect('equal', 'box') 

    ax1.set_title('(a) PVCG Payment', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = -0.5)
    ax2.set_title('(b) Procurement Level', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = -0.5)

    ax1.set_xlim(0., 250.)
    ax1.set_ylim(0., 250.)
    ax1.set_xticks(np.arange(0., 251., 50.))
    ax1.set_yticks(np.arange(0., 251., 50.))

    ax2.set_xlim(0., 5.5)
    ax2.set_ylim(0., 5.5)
    ax2.set_xticks(np.arange(0., 6., 1.))
    ax2.set_yticks(np.arange(0., 6., 1.))

    ax1.plot(x_grid, Y_grid, color = 'grey', linestyle = ':', linewidth = 0.75)
    ax2.plot(x_grid, Y_grid, color = 'grey', linestyle = ':', linewidth = 0.75)

    for i in range(N):
        ax1.scatter(p_vec_closed[i], payments_rim[i], 
                    color = COLOR_PALETTE[i], marker = MARKER_STYLE,
                    s = MARKER_SIZE_BIG,
                    label = 'Supplier ' + str(i+1))
        ax2.scatter(x_vec_closed[i], x_vec_rim[i], 
                    color = COLOR_PALETTE[i], marker = MARKER_STYLE,
                    s = MARKER_SIZE_BIG)
    
    ax3.set_aspect('equal', 'box')
    # print(z_diff)
    sns.heatmap(z_diff, cmap=my_cmap, ax = ax3,
                annot=False, linewidths=0.5, center = 0.,
                vmin = -0.035, vmax = 0.035,
                fmt = '.0%',
                cbar_kws={
                    "orientation": "vertical",
                    "location": "right" ,
                    "shrink": 0.52,
                    "format": PercentFormatter(1, decimals=0)}
                )
    cbar = ax3.collections[0].colorbar
    cbar.set_ticks([-0.03, -0.02,  -0.01, 
                    0., 0.01, 0.02, 0.03])

    ax3.set_xticklabels(['1', '2', '3', '4', '5',
                         '6', '7', '8', '9', '10'])
    ax3.set_yticklabels(['1', '2', '3', '4', '5',
                         '6', '7', '8', '9', '10'])
    ax3.invert_yaxis()
    ax3.set_xlabel('Supplier ID', size = LABEL_SIZE_SINGLE_COLUMN,
                    labelpad= 0)
    ax3.set_ylabel('Supplier ID', size = LABEL_SIZE_SINGLE_COLUMN)   
    ax3.set_title('         (c) Relative Error of z-Matrix', 
                size = TITLE_SIZE_SINGLE_COLUMN,
                y = -0.46)

    # fig.tight_layout()
    fig.subplots_adjust(top = 1., hspace=0.5, wspace=0.25)

    fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.), ncol=5)

    plt.savefig(PIC_NAME_SCATTER, 
                transparent=True, bbox_inches='tight', pad_inches=0.02)

    # plt.show()


if __name__ == '__main__':

    plot_scatter()