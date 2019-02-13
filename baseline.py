import get_data
import os
import numpy as np
import matplotlib.pyplot as plt


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    
    xsamp=box_pts*np.arange(int(len(y)/box_pts))
    return y_smooth[xsamp]

# ---------------------------------------------------------------------------
def check_baseline(fulldirname, config):

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_BASELINE")

    i_test_deb, i_test_fin = 21, 40
    test = "IQ-ALL_Science-Data"
    fichlist = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat' \
                and f[i_test_deb:i_test_fin]==test]

    # -----------------------------------------------------------------------
    # Processing baseline 
    if len(fichlist)>0:
        Chan0_i, Chan0_q, _, _, _ = get_data.readIQ(os.path.join(datadirname,fichlist[0]))
        l=len(Chan0_i[:,0])
        mod=np.zeros((l,41))
        for pix in range(41):
            mod[:,pix]=np.sqrt(Chan0_i[:,pix].astype('float')**2 + Chan0_q[:,pix].astype('float')**2)

        Chan0_i, Chan0_q=0,0

        t = np.arange(l)/(20e6/2**7)
        n_boxes=41
        n_lines=6
        n_cols=7

        fig = plt.figure(figsize=(18, 12))
        for box in range(n_boxes):
            marge = 0.5
            ymax = mod[:,box].max() + (mod[:,box].max() - mod[:,box].min())*marge
            ymin = mod[:,box].min() - (mod[:,box].max() - mod[:,box].min())*marge
            ax = fig.add_subplot(n_lines, n_cols, box+1)
            ax.plot(t, mod[:,box])
            ax.set_ylim([ymin, ymax])
            ax.set_title(r'Pixel {0:2d}'.format(box))
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            ratio = 100/2**15
            ax2 = plt.gca().twinx()
            ax2.set_ylim([ymin*ratio, ymax*ratio])
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            if box//n_cols == n_lines-1:
                ax.set_xlabel(r'Time (s)')
            else:
                plt.xticks(visible=False)
            if box%n_cols == 0:
                ax.set_ylabel(r'Module (A.U.)')
            if box%n_cols == n_cols-1:
                ax2.set_ylabel(r'Module (% of FSR)')

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ax2.yaxis.label]):
                item.set_fontsize(10)
                item.set_weight('bold')
            for item in (ax.get_xticklabels() + ax.get_yticklabels() + ax2.get_yticklabels()):
                item.set_fontsize(8)
            for item in (ax.get_yticklabels() + ax2.get_yticklabels()):
                item.set_rotation(45)

        fig.tight_layout()
        plt.savefig(pltfilename+'_40pix.png', bbox_inches='tight')

# ---------------------------------------------------------------------------
