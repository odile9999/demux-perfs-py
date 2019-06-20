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
def meas_energy_r(fulldirname, config, pix=40, export_npy_files=False):

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

        time_sec = time_stamps / (config['fs']/2**config['power_to_fs2'])

        # exporting data in npy file
        if export_npy_files:
            np.save(os.path.join(datadirname,'events.npy'),([time_sec, baseline, energy]))

        # Removing huge errors if any
        i_bad = np.where(abs(energy-energy.mean()) > 10*energy.std())[0]
        for i in i_bad:
            print("WARNING !! removing energy at index {0:6d}. Probably a wrong measurment.".format(i))
            energy=np.delete(energy, i)
            time_sec=np.delete(time_sec, i)
            baseline=np.delete(baseline, i)
        if len(i_bad)>0:
            pltfilename = pltfilename+"_ERROR_"


        # Making the histogram plot
        nrj, nrj_resol_at_7kev = fit_tools.gauss_fit(energy,fit_tools.number_of_bins(energy),show=True, \
            pltfilename=pltfilename+'_'+str(file_counter), inf=None)
        tab_nb_events = np.append(tab_nb_events, len(time_sec))
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

# -----------------------------------------------------------------------------
def get_pix_iq(filename, config, pix):
    fs=config["fs"]/2**config["power_to_fs2"]
    print(filename)
    chan0_i, chan0_q, _, _, _ = get_data.read_iq(filename)
    module = np.sqrt(chan0_i[:, pix].astype(float)**2 + chan0_q[:, pix].astype(float)**2)
    t=np.arange(len(module))/fs
    return(t, module)

# -----------------------------------------------------------------------------
def make_npy(fulldirname, config, pix):

    dirnm = os.path.join(fulldirname, config['dir_data'])

    noise_file = [f for f in os.listdir(dirnm) if os.path.isfile(os.path.join(dirnm,f)) and f[-18:]=='-Noise-Charact.dat']
    pulse_file = [f for f in os.listdir(dirnm) if os.path.isfile(os.path.join(dirnm,f)) and f[-23:]=='-Pulse-Charact@7keV.dat']
    meas_file = [f for f in os.listdir(dirnm) if os.path.isfile(os.path.join(dirnm,f)) and f[:9]=='dump_dre_']

    if len(noise_file)>0 and len(pulse_file)>0 and len(meas_file):
        noise_t, noise_a = get_pix_iq(os.path.join(dirnm,noise_file[0]), config, pix)
        pulse_t, pulse_a = get_pix_iq(os.path.join(dirnm,pulse_file[0]), config, pix)
        measure_t, measure_a = get_pix_iq(os.path.join(dirnm,meas_file[0]), config, pix)

        np.save(os.path.join(dirnm,'noise.npy'),([noise_t, noise_a]))
        np.save(os.path.join(dirnm,'pulse.npy'),([pulse_t, pulse_a]))
        np.save(os.path.join(dirnm,'measure.npy'),([measure_t, measure_a]))

    else:
        print("Missing file for the export of EP data.")

# -----------------------------------------------------------------------------

def check_npy_files(fulldirname, config):

    dirnm = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    general_tools.checkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "CHECK_NPY_FILES")
    
    events_file = os.path.join(dirnm,'events.npy')
    noise_file = os.path.join(dirnm,'noise.npy')
    pulse_file = os.path.join(dirnm,'pulse.npy')
    measure_file = os.path.join(dirnm,'measure.npy')

    if os.path.isfile(events_file) and os.path.isfile(noise_file) and \
        os.path.isfile(pulse_file) and os.path.isfile(measure_file):
        t_event, baseline, energy = np.load(os.path.join(dirnm,'events.npy'))
        noise_t, noise_a = np.load(os.path.join(dirnm,'noise.npy'))
        pulse_t, pulse_a = np.load(os.path.join(dirnm,'pulse.npy'))
        measure_t, measure_a = np.load(os.path.join(dirnm,'measure.npy'))

        # Making the plot
        fig = plt.figure(figsize=(7, 14))

        ax1 = fig.add_subplot(5, 1, 1)
        ax1.plot(t_event, baseline)
        ax1.set_title('Baseline')
        ax1.set_xlabel('time (sec)')

        ax2 = fig.add_subplot(5, 1, 2)
        ax2.plot(t_event, energy)
        ax2.set_title('Energies')
        ax2.set_xlabel('time (sec)')

        ax3 = fig.add_subplot(5, 1, 3)
        ax3.plot(noise_t, noise_a)
        ax3.set_title('Noise data')
        ax3.set_xlabel('time (sec)')

        ax4 = fig.add_subplot(5, 1, 4)
        ax4.plot(pulse_t, pulse_a)
        ax4.set_title('Pulse data')
        ax4.set_xlabel('time (sec)')

        ax5 = fig.add_subplot(5, 1, 5)
        ax5.plot(measure_t, measure_a)
        ax5.set_title('Measurement data')
        ax5.set_xlabel('time (sec)')

        fig.tight_layout()
        plt.savefig(pltfilename+'.png', bbox_inches='tight')

