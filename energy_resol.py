# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:28:01 2018

@author: Laurent
"""

import numpy as np
import os
import get_data, fit
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------
def meas_energy_r(fulldirname, config):

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_ENERGY-RESOL")

    # Reading data from files    
    events_name = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f=='events.dat']
    
    if len(events_name)>0:
        time_stamps, _, energy, baseline = get_data.readEvents(os.path.join(datadirname, events_name[0]))
        time_sec = time_stamps / (config['fs']/config['power_to_fs2'])
        # Making the histogram plot
        fit.hist_and_fit(energy,fit.number_of_bins(energy),show=True, pltfilename=pltfilename, inf=None, out=True)

        # Making the baseline plot
        ymax = np.max(baseline) + (np.max(baseline) - np.min(baseline))*0.2
        ymin = np.min(baseline) - (np.max(baseline) - np.min(baseline))*0.2
        fig = plt.figure(figsize=(9, 5))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(time_sec, baseline)
        ax.set_xlabel(r'Time (s)')
        ax.set_ylabel(r'Baseline (A.U.)')
        ratio = 100/2**15
        ax2 = plt.gca().twinx()
        ax2.set_ylim([ymin*ratio, ymax*ratio])
        ax2.set_ylabel(r'Baseline (% of FSR)')
        for item in (ax.yaxis.label, ax.xaxis.label, ax2.yaxis.label):
            item.set_weight('bold')
            item.set_fontsize(15)

        fig.tight_layout()
        plt.savefig(pltfilename+'_BASELINE.png', bbox_inches='tight')
   

#------------------------------------------------------------------------------
def plot_GSE_spectrum(fulldirname, config, Emin, Emax):

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_GSE_SPECTRUM")

    # Reading data from files
    i_type_fin = 9
    spectrum_name = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[:i_type_fin]=='Spectrum_']
    hysto_en, hysto_counts = get_data.readSpectrum(os.path.join(datadirname, spectrum_name[0]))
    
    # Making the plot
    fig = plt.figure(figsize=(9, 5))
    MaxCounts = np.max(hysto_counts)*1.3

    ax7 = fig.add_subplot(1, 1, 1)
    ax7.plot(hysto_en, hysto_counts, 'o')
    ax7.set_ylabel(r'Counts')
    ax7.set_xlabel(r'Energy (eV)')
    ax7.axis([Emin, Emax, 0, MaxCounts])

    fig.tight_layout()
    plt.savefig(pltfilename+'.png', bbox_inches='tight')

    #------------------------------------------------------------------------------
