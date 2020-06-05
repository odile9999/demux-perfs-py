# -----------------------------------------------------------------------------
# Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import get_data, general_tools, fit_tools

# -----------------------------------------------------------------------
def get_files_freq(fulldirname, dumptype):
    r"""
        This function extract file names and frequency values from 
        a directory containing data files.
        
        Parameters
        ----------
        dirname : string
        The name of the directory containing the data files

        dumptype : defines the kind of files we are interested in.
        "IQ-ALL" or "IQ-TES"

        Returns
        ------- 
        fichlist: array
        Contains the file names

        freqs: array
        Contains the frequencies found in the file names

        """
    i_fich_deb, i_fich_fin = 16, 20
    i_type_deb, i_type_fin = 21, 27
    i_freq_deb, i_freq_fin = 32, 39
    i_gbw_deb, i_gbw_fin = 27, 31
    
    fnames = [f for f in os.listdir(fulldirname) \
                if os.path.isfile(os.path.join(fulldirname, f)) \
                and f[-4:]=='.dat'\
                and f[i_type_deb:i_type_fin]==dumptype\
                and f[i_gbw_deb:i_gbw_fin]=="_GBW"]
 
    if len(fnames)>0:
        freqs = file_nb = np.zeros(len(fnames))
        fichlist = fnames.copy()
        for i in range(len(fnames)):
            name=fnames[i]
            file_nb=int(name[i_fich_deb: i_fich_fin])
            fichlist[file_nb]=name 
            freqs[file_nb]=float(name[i_freq_deb:i_freq_fin])
    else:
        fichlist=freqs=[] 
    return(fichlist, freqs)
        
# -----------------------------------------------------------------------
def get_bw(amplitudes_db, level_db=3):
    r"""
        This function measures the bandwith of a low pass band at a
        specified level.
        
        Parameters
        ----------
        amplitudes_db : array
        contains the values of the amplitudes.

        level_db : number
        level in dB to define the out of band limit (default = 3dB)

        Returns
        ------- 
        i_max : number
        the index of the last "in-band amplitude"
        """
    
    i_inband = np.where(amplitudes_db >= np.max(amplitudes_db-1*level_db))
    i_max = i_inband[0][-1]
    return(i_max)

# -----------------------------------------------------------------------
def color_ok(val, val_min, val_max):
    r"""
        This function checks a value wrt to levels and returns a color
        indicator: 'g' (green) if the value is between the two levels
        'r' (red) if not.
    """
    if val > val_min and val < val_max:
        color_ok='g'
    else:
        color_ok = 'r'
    return(color_ok)

# -----------------------------------------------------------------------
def get_gbw_freq_and_amp(datadirname, dumptype, chan):
    r"""
        This function reads data from several DRE data files and measures
        the amplitude of the signal as a function of the frequency.
        
        Parameters
        ----------
        datadirname : string
        The name of the directory containing the data files

        dumptype : defines the kind of files we are interested in.
        "IQ_ALL" or "IQ_TES"

        chan : number
        Specifies the channel to be processed (0 or 1)

        Returns
        ------- 
        frequency: array
        frequencies used for the test

        amplitude:
        amplitude measured at each frequency 

        gain_dre:
        value of the DRE gain parameter
        """
    pix_test=40
    fichlist, f = get_files_freq(datadirname, dumptype)
    i_gain_deb, i_gain_fin = 46, -4

    amplitudes, gain_dre = [], 0

    if len(fichlist) !=0:
        gain_dre = fichlist[0][i_gain_deb: i_gain_fin]
        nfiles = len(fichlist)
        print("{0:3d} files to process...".format(nfiles))

        for fich in fichlist:
            print('Processing file ' + fich)
            if dumptype=="IQ-ALL":
                chan0_i, chan0_q, chan1_i, chan1_q, _ = get_data.read_iq(os.path.join(datadirname, fich))
                if chan==0:
                    modulus = \
                        np.sqrt(chan0_i[:,pix_test].astype('float')**2 + chan0_q[:,pix_test].astype('float')**2)
                else:
                    modulus = \
                        np.sqrt(chan1_i[:,pix_test].astype('float')**2 + chan1_q[:,pix_test].astype('float')**2)    
            else: # Le dumptype est "IQ-TEST"
                data, _ = get_data.readfile(os.path.join(datadirname, fich))
                modulus = np.sqrt(data[1:,0].astype('float')**2 + data[1:,1].astype('float')**2)

            amplitudes=np.append(amplitudes, modulus.max()-modulus.min())
    
    return(f, amplitudes, gain_dre)
                    
# -----------------------------------------------------------------------
def process_gbw(fulldirname, config, chan, margin_pc=5):
    r"""
        This function reads data from several DRE IQ data files
        to compute the Gain BadWidth product (GBW).
        
        Parameters
        ----------
        fulldirname: string
        The name of the directory containing the data files

        config: dictionnary
        Contains path and constants definitions

        chan: number
        Specifies the channel to be processed (0 or 1)
    
        margin_pc: number
        Acceptable margin, as a percentage, on the value of the GBWP (default is 5) 

        Returns
        ------- 
        gbwp_ok: boolean
        True if gain bandwidth product is in the correct range.
        
        """
    fs = config["fs"]

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    general_tools.checkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_GBW.png")

    f, a, gain_dre = get_gbw_freq_and_amp(datadirname, "IQ-TST", chan)

    if len(f) > 0:
        a_db = 20*np.log10(a/a.max())

        # Fit of the transfer function and measurement of the GBP
        fit_params=fit_tools.low_pass_fit(f, a)
        f_high_res = np.linspace(1, 1e6, num=100000)
        fit = fit_tools.low_pass(f_high_res, fit_params[0], fit_params[1])
        fit_db = 20*np.log10(fit/fit.max())
        gbw = 1./fit_params[1]
        gbwp_ok = (gbw>=config['GBWP']*(1-margin_pc/100)) and (gbw<=config['GBWP']*(1+margin_pc/100))

        min_amp_db, max_amp_db = -80, 10
        gbw_min, gbw_max = 13e3, 17e3
        fs_sur2 = fs/(2*2**7)

        gbw_color = color_ok(gbw, gbw_min, gbw_max)

        # Plot
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.semilogx(f, a_db, 's', label='Transfer function before decimation filter')
        ax.semilogx(f_high_res, fit_db, color='k', label='Fit (LPF 1st order)')

        # Processing filtered data
        f2, a2, _ = get_gbw_freq_and_amp(datadirname, "IQ-ALL", chan)
        if len(f2) > 0:
            a2_db = 20*np.log10(a2/a2.max())
            ax.semilogx(f2, a2_db, '.', label='Transfer function after decimation filter')

        ax.semilogx([gbw, gbw], [min_amp_db, max_amp_db], '--', color='b', label='BBFB Gain Bandwidth Product')

        # Adding some references and text on the plot
        comment1 = 'Channel {0:1d}'.format(chan)
        comment2 = 'GainDRE = ' + gain_dre
        comment3 = 'GBW = {0:6.0f} Hz'.format(gbw)
        ax.semilogx([fs_sur2, fs_sur2], [min_amp_db, max_amp_db], '--', color='k' , linewidth=0.5 , label='Fs / 2')
        ax.text(20, min_amp_db + 35, comment1, color='b')
        ax.text(20, min_amp_db + 30, comment2, color='b')
        ax.text(20, min_amp_db + 25, comment3, color=gbw_color)

        ax.legend()
        ax.axis([10, 1e6, min_amp_db, max_amp_db])
        ax.set_xlabel(r'Frequency (Hz)')
        ax.set_ylabel(r'Amplitude (dB)')
        ax.set_title(r'Bandwidth measurement')
        major_ticks = np.arange(min_amp_db, max_amp_db, 10)
        minor_ticks = np.arange(min_amp_db, max_amp_db, 2)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        fig.tight_layout()
        plt.savefig(pltfilename, bbox_inches='tight')
    else:
        f=a_db=0
        gbwp_ok=False

    return(gbwp_ok)

# -----------------------------------------------------------------------
