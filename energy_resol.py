# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:28:01 2018

@author: Laurent
"""

import numpy as np
import os
import get_data, fit_tools, general_tools
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------
def meas_energy_r(fulldirname, config, pix=40):

    session_info = general_tools.get_session_info(fulldirname)
    if 'BackupVersion' in session_info:
        backup_version = np.float(session_info['BackupVersion'])
    else:
        backup_version = 0

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    general_tools.checkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_ENERGY-RESOL")
    logfilename = os.path.join(plotdirname, "ENERGY-RESOL.log")

    # Reading data from files
    if backup_version < 1:
        events_name = [f for f in os.listdir(datadirname) \
                    if os.path.isfile(os.path.join(datadirname, f)) \
                    and os.path.getsize(os.path.join(datadirname, f))!=0 \
                    and f[:6]=='events' and f[-4:]=='.dat']
    else:
        events_name = [f for f in os.listdir(datadirname) \
                    if os.path.isfile(os.path.join(datadirname, f)) \
                    and os.path.getsize(os.path.join(datadirname, f))!=0 \
                    and f[-7:]=='.events']

    tab_nb_events = np.array([], dtype=int)
    tab_nrj = tab_nrj_resol_at_7kev = np.array([])
    file_counter = 0
    for file in events_name:
        time_stamps, _, pix_id, energy, baseline = get_data.read_events(os.path.join(datadirname, file), backup_version)

        # keeping only pixels from the pixel we are interested in
        i_good = np.where(pix_id==pix)
        time_stamps=time_stamps[i_good[0]]
        energy=energy[i_good[0]]
        baseline=baseline[i_good[0]]

        # Removing huge errors if any
        i_bad = np.where(abs(energy-energy.mean()) > 10*energy.std())[0]
        for i in i_bad:
            print("WARNING !! removing energy at index {0:6d}. Probably a wrong measurment.".format(i))
            energy=np.delete(energy, i)
            time_stamps=np.delete(time_stamps, i)
            baseline=np.delete(baseline, i)
        if len(i_bad)>0:
            pltfilename = pltfilename+"_ERROR_"

        time_sec = time_stamps / (config['fs']/2**config['power_to_fs2'])

        # Making the histogram plot
        nrj, nrj_resol_at_7kev = fit_tools.gauss_fit(energy,fit_tools.number_of_bins(energy),show=True, \
            pltfilename=pltfilename+'_'+str(file_counter), inf=None)
        tab_nb_events = np.append(tab_nb_events, len(time_stamps))
        tab_nrj = np.append(tab_nrj, nrj)
        tab_nrj_resol_at_7kev = np.append(tab_nrj_resol_at_7kev, nrj_resol_at_7kev)

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
        plt.savefig(pltfilename+'_'+str(file_counter)+'_BASELINE.png', bbox_inches='tight')
        file_counter+=1

    if file_counter>0:
        flog = open(logfilename, 'w')
        flog.write("/===============================================================|\n")
        flog.write("Results of energy resolution measurements:\n")
        flog.write("Number of event lists processed: {0:4d}\n".format(file_counter))
        flog.write("Number of events:\n")
        flog.write("   Minimum: {0:6d}\n".format(tab_nb_events.min()))
        flog.write("   Maximum: {0:6d}\n".format(tab_nb_events.max()))
        flog.write("   Mean:      {0:6.3f}\n".format(tab_nb_events.mean()))
        flog.write("Mean of the mean energy:             {0:8.3f} eV\n".format(tab_nrj.mean()))
        flog.write("Mean of the energy resolution:       {0:8.3f} eV\n".format(tab_nrj_resol_at_7kev.mean()))
        if file_counter>1:
            flog.write("Dispersion of the mean energy:       {0:8.3f} eV\n".format(tab_nrj.std()))
            flog.write("Dispersion of the energy resolution: {0:8.3f} eV\n".format(tab_nrj_resol_at_7kev.std()))
        flog.write("/===============================================================|\n")
        flog.close()

    return(tab_nb_events, tab_nrj, tab_nrj_resol_at_7kev)

#------------------------------------------------------------------------------
def plot_gse_spectrum(fulldirname, config, e_min, e_max):

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    general_tools.checkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_GSE_SPECTRUM")

    # Reading data from files
    i_type_fin = 9
    spectrum_name = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[:i_type_fin]=='Spectrum_']
    hysto_en, hysto_counts = get_data.read_spectrum(os.path.join(datadirname, spectrum_name[0]))

    # Making the plot
    fig = plt.figure(figsize=(9, 5))
    max_counts = np.max(hysto_counts)*1.3

    ax7 = fig.add_subplot(1, 1, 1)
    ax7.plot(hysto_en, hysto_counts, 'o')
    ax7.set_ylabel(r'Counts')
    ax7.set_xlabel(r'Energy (eV)')
    ax7.axis([e_min, e_max, 0, max_counts])

    fig.tight_layout()
    plt.savefig(pltfilename+'.png', bbox_inches='tight')

    #------------------------------------------------------------------------------
