"""Plot the Fig. 1 in the paper.
"""
import csv
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns


PIC_NAME = 'fig1_double.png'
INPUT_NAME = 'metric_csv.csv'

FIG_SIZE_SINGLE_COLUMN = (8, 4)
LABEL_SIZE_SINGLE_COLUMN = 10
# TITLE_SIZE_SINGLE_COLUMN = 16
TITLE_SIZE_SINGLE_COLUMN = 12

# FIG_SIZE_DOUBLE_COLUMN = (16, 4)
# LABEL_SIZE_DOUBLE_COLUMN = 10
# TITLE_SIZE_DOUBLE_COLUMN = 12

LABEL_SIZE_3D = 10
TITLE_SIZE_3D = 12

LINE_STYLE = '-'
LINE_WIDTH = 1.


COLOR_PALETTE = ["#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b",
                 "#bdcf32", "#87bc45", "#00bfa0", "#27aeef", "#6481D9"]
MARKER_STYLE = "s"
MARKER_SIZE = 4
MARKER_SIZE_BIG = 30

Y_SUBTITLE = -0.4

COLOR = 'blue'
LINE_STYLE = ':'
LINE_WIDTH = 0.75

my_cmap = sns.diverging_palette(120., 255., s=100, as_cmap=True)

def plot_scatter():
    """The scatter plot for the RIM results and the closed-form solutions.
    """

    n_label = []
    capital_n_label = []
    n_series = []
    capital_n_series = []
    setup_time = []
    srs_size = []
    commit_time = []
    prover_time = []
    transcript_size = []
    verifier_time = []

   
    with open(INPUT_NAME, 'r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)

        for row in reader:

            n_label.append(row[0])
            capital_n_label.append(row[1])
            n_series.append(int(row[2]))
            capital_n_series.append(int(row[3]))
            setup_time.append(float(row[4]))
            srs_size.append(float(row[5]))
            commit_time.append(float(row[6]))
            prover_time.append(float(row[7]))
            transcript_size.append(float(row[8])/1000.)
            verifier_time.append(float(row[9]))

    x_tick_label = [n_label[i] for i in range(0, len(n_label), 2)]
    x_tick = [n_series[i] for i in range(0, len(n_series), 2)]
    
    x_capital_tick_label = [capital_n_label[i] for i in range(0, len(n_label), 2)]
    x_capital_tick = [capital_n_series[i] for i in range(0, len(n_series), 2)]

    fig = plt.figure(figsize=(6,9))
    gs = gridspec.GridSpec(
        3, 2, 
        width_ratios=(1.,1.), wspace=0.5, hspace=0.6) 

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = plt.subplot(gs[3])
    ax5 = plt.subplot(gs[4])
    ax6 = plt.subplot(gs[5])

    # ax1.set_aspect('equal', 'box') 
    # ax2.set_aspect('equal', 'box')
    # ax3.set_aspect('equal', 'box') 
    # ax4.set_aspect('equal', 'box')
    # ax5.set_aspect('equal', 'box') 
    # ax6.set_aspect('equal', 'box')

    ax1.set_title('(a) Setup Time (log-log)', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = 1.1)
    ax2.set_title('(b) SRS Size (log-log)', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = 1.1)
    ax3.set_title('(c) Commitment Time (log-log)', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = 1.1)
    ax4.set_title('(d) Prover Time (log-log)', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = 1.1)
    ax5.set_title('(e) Transcript Size (linear-log)', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = 1.1)
    ax6.set_title('(f) Verifier Time (linear-log)', size = TITLE_SIZE_SINGLE_COLUMN,
                  y = 1.1)


    ax1.set_xlabel('Matrix Dimension (n)', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    ax1.set_ylabel('Time (seconds)', size = LABEL_SIZE_SINGLE_COLUMN)  
    ax1.loglog(n_series, setup_time, 
               marker = MARKER_STYLE, markersize = MARKER_SIZE,
                 color = COLOR, alpha = 0.75, 
                 linestyle = LINE_STYLE, linewidth = LINE_WIDTH)
    ax1.grid(True, which="major", ls="-", color='0.93')
    ax1.set_xticks(x_tick, x_tick_label)

    ax2.set_xlabel('Matrix Dimension (n)', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    ax2.set_ylabel('Size (bytes)', size = LABEL_SIZE_SINGLE_COLUMN)  
    ax2.loglog(n_series, srs_size, 
               marker = MARKER_STYLE, markersize = MARKER_SIZE,
                 color = COLOR, alpha = 0.75, 
                 linestyle = LINE_STYLE, linewidth = LINE_WIDTH)
    ax2.grid(True, which="major", ls="-", color='0.93')
    ax2.set_xticks(x_tick, x_tick_label)

    ax3.set_xlabel('# of Non-Zero Entries (N/3)', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    ax3.set_ylabel('Time (seconds)', size = LABEL_SIZE_SINGLE_COLUMN)  
    ax3.loglog(capital_n_series, srs_size, 
               marker = MARKER_STYLE, markersize = MARKER_SIZE,
                 color = COLOR, alpha = 0.75, 
                 linestyle = LINE_STYLE, linewidth = LINE_WIDTH)
    ax3.grid(True, which="major", ls="-", color='0.93')
    ax3.set_xticks(x_capital_tick, x_capital_tick_label)

    ax4.set_xlabel('# of Non-Zero Entries (N/3)', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    ax4.set_ylabel('Time (seconds)', size = LABEL_SIZE_SINGLE_COLUMN)  
    ax4.loglog(capital_n_series, prover_time, 
               marker = MARKER_STYLE, markersize = MARKER_SIZE,
                 color = COLOR, alpha = 0.75, 
                 linestyle = LINE_STYLE, linewidth = LINE_WIDTH)
    ax4.grid(True, which="major", ls="-", color='0.93')
    ax4.set_xticks(x_capital_tick, x_capital_tick_label)
    
    ax5.set_xlabel('Matrix Dimension (n)', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    ax5.set_ylabel('Size (KB)', size = LABEL_SIZE_SINGLE_COLUMN)  
    ax5.semilogx(n_series, transcript_size, 
                 marker = MARKER_STYLE, markersize = MARKER_SIZE,
                 color = COLOR, alpha = 0.75, 
                 linestyle = LINE_STYLE, linewidth = LINE_WIDTH)
    ax5.grid(True, which="major", ls="-", color='0.93')
    ax5.set_xticks(x_tick, x_tick_label)


    ax6.set_xlabel('Matrix Dimension (n)', size = LABEL_SIZE_SINGLE_COLUMN,
                   labelpad=4)
    ax6.set_ylabel('Time (milliseconds)', size = LABEL_SIZE_SINGLE_COLUMN)  
    ax6.semilogx(n_series, verifier_time, 
                 marker = MARKER_STYLE, markersize = MARKER_SIZE,
                 color = COLOR, alpha = 0.75, 
                 linestyle = LINE_STYLE, linewidth = LINE_WIDTH)
    ax6.grid(True, which="major", ls="-", color='0.93')
    ax6.set_xticks(x_tick, x_tick_label)

    # fig.tight_layout()
    # fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.), ncol=5)

    plt.savefig(PIC_NAME, 
                transparent=True, bbox_inches='tight', pad_inches=0.02)

    plt.show()


if __name__ == '__main__':

    plot_scatter()
