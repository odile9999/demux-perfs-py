# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
"""
    scan_feedback_tools module
    ====================

    Developped by: L. Ravera

    Project: Athena X-IFU / DRE-DEMUX

    General purpose tools needed to read and plot the scan feedback data

    """

# -----------------------------------------------------------------------
# Imports
import os, csv
import general_tools
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------
def get_scanfb(fulldirname, config):
    r"""
        This function gets the scan feedback data from a csv file.

        Parameters:
        -----------
        fulldirname: string
        The name of the scanfeedback data file (with the path)

        config: dictionnary
        Contains path and constants definitions

        Returns
        -------
        scanfbdata : dictionnary

        """
    
    scanfbdata={}
    dirname = os.path.join(os.path.normcase(fulldirname), os.path.normcase(config['dir_data']))
    scanfbfullfilename=os.path.join(dirname, 'scanFB.dat')
    if os.path.isfile(scanfbfullfilename):
        print("Reading Scan-feedback data...")
        with open(scanfbfullfilename, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=';', quotechar='|')

            for row in reader:
                n = len(row)
                if reader.line_num == 1:
                    keys = row
                    for i in range(n):
                        scanfbdata[keys[i]]=np.array([]) # Initialises an empty array
                else:
                    for i in range(n):
                        if keys[i] == 'Date':
                            scanfbdata[keys[i]] = np.append(scanfbdata[keys[i]], row[i])
                        else:
                            scanfbdata[keys[i]] = np.append(scanfbdata[keys[i]], float(row[i].replace(',','.')))
    else:
        print("Scan-feedback file not found.")
        scanfbdata=0

    return(scanfbdata)

# -----------------------------------------------------------------------
def plot_scanfb(scanfbdata, fulldirname, config, limit_error):
    r"""
        This function plots scan feedback data.

        Parameters
        ----------
        scanfbdata: Dictionnary
        scan feedback data (frequencies, module, phase minus best fit) values.

        fulldirname: string
        The name of the dump file (with the path)

        config: dictionnary
        Contains path and constants definitions

        limit_error: number
        Acceptance level for the phase-fit values (default=1.5 degrees)

        Returns
        -------
        Nothing

        """

    fmin_khz=1e3
    fmax_khz=5e3

    plotdirname = os.path.join(os.path.normcase(fulldirname), os.path.normcase(config['dir_plots']))
    general_tools.checkdir(plotdirname)
    pltfilename = os.path.join(plotdirname, 'PLOT_SCANFEEDBACK.png')

    fig = plt.figure(figsize=(12, 8))

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(scanfbdata['freq(kHz)'], scanfbdata['Mod(dB)'], linewidth=3)
    axes=plt.gca()
    ymin, ymax=axes.get_ylim()
    ax1.plot([fmin_khz, fmin_khz], [ymin, ymax], '--r', linewidth=1)
    ax1.plot([fmax_khz, fmax_khz], [ymin, ymax], '--r', linewidth=1)
    ax1.set_title('Scan-Feedback plot')
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    #ax1.set_xlabel('Frequency (kHz)')
    ax1.set_xticklabels([])
    ax1.set_ylabel('Module (dB)')

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(scanfbdata['freq(kHz)'], scanfbdata['Phi-fit(deg)'], linewidth=3)
    axes=plt.gca()
    ymin, ymax=axes.get_ylim()
    ax2.plot([fmin_khz, fmin_khz], [-1*limit_error, limit_error], '--r', linewidth=1)
    ax2.plot([fmax_khz, fmax_khz], [-1*limit_error, limit_error], '--r', linewidth=1)
    ax2.plot([fmin_khz, fmax_khz], [limit_error, limit_error], '--r', linewidth=1)
    ax2.plot([fmin_khz, fmax_khz], [-1*limit_error, -1*limit_error], '--r', linewidth=1)
    xmin, xmax=scanfbdata['freq(kHz)'].min(), scanfbdata['freq(kHz)'].max()
    ax2.plot([xmin, xmax], [0, 0], '-k', linewidth=2)
    ax2.grid(color='k', linestyle=':', linewidth=0.5)
    ax2.set_xlabel('Frequency (kHz)')
    ax2.set_ylabel('Phase - fit (deg)')
    fig.tight_layout()
    plt.savefig(pltfilename, bbox_inches='tight')

# -----------------------------------------------------------------------
def check_scanfb(fulldirname, config, limit_error=1.5, make_plot=True):
    r"""
        This function reads, checks and plots the data of a scanfeedback.

        Parameters
        ----------
        fulldirname: string
        The name of the scanfeedback data file (with the path)

        config: dictionnary
        Contains path and constants definitions

        limit_error: number
        Acceptance level for the phase-fit values (default=1.5 degrees)

        make_plot: boolean
        If True the function makes a plot with the feedback data.

        Returns
        -------
        Nothing

        scanfb_ok: boolean
        True if scanfeedback is considered as successfull
        """
    fmin_khz=1e3
    fmax_khz=5e3

    scanfb_ok=False
    sfb_dat = get_scanfb(fulldirname, config)

    if sfb_dat !=0:
        if make_plot:
            plot_scanfb(sfb_dat, fulldirname, config, limit_error)
        # computing index of the band of interest
        i_interest=np.where((sfb_dat['freq(kHz)']>fmin_khz) & (sfb_dat['freq(kHz)']<fmax_khz))[0]
        scanfb_ok=abs(sfb_dat['Phi-fit(deg)'][i_interest]).max()<=limit_error
        if not scanfb_ok:
            print('   >> Warning Phase compensation error is greater than {0:3.1} deg'.format(limit_error))
    else:
        print('  >> Warning! There is no scanfeedback data.')

    return(scanfb_ok)

# -----------------------------------------------------------------------
