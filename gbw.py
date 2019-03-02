# -----------------------------------------------------------------------------
# Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import get_data, general_tools

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
        fichlist=freqs=0 
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
def process_gbw(fulldirname, config, dumptype, chan):
    r"""
        This function reads data from several DRE IQ data files
        to compute the Gain BadWidth product (GBW).
        
        Parameters
        ----------
        fulldirname : string
        The name of the directory containing the data files
        (no path)

        config : dictionnary
        Contains path and constants definitions

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
        """
    fs = config["fs"]
    pix_test=40

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    general_tools.checkdir(plotdirname)

    if dumptype == "IQ-ALL":
        pltfilename = os.path.join(plotdirname, "PLOT_GBW_FILTERED.png")
    else:
        pltfilename = os.path.join(plotdirname, "PLOT_GBW_UNFILTERED.png")

    fichlist, f = get_files_freq(datadirname, dumptype)

    if fichlist !=0:
        nfiles = len(fichlist)
        amplitudes = np.zeros(nfiles)

        print("{0:3d} files to process...".format(nfiles))
        file_nb=0

        for fich in fichlist:
            print('Processing file {0:3d}/{1:3d} '.format(file_nb, nfiles), end='')
            print(fich)
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

            amplitudes[file_nb]=np.max(modulus)-np.min(modulus)
                    
            file_nb+=1

        a_db = 20*np.log10(amplitudes/np.max(amplitudes))
        # Measuring 3dB BW
        i_bw = get_bw(a_db, 3)
        gbw = f[i_bw]
        gbw_min = 12e3
        gbw_max = 18e3
        fs_sur2 = fs/(2*2**7)

        i_gain_deb, i_gain_fin = 46, -4
        gain_dre = fichlist[0][i_gain_deb: i_gain_fin]

        gbw_color = color_ok(gbw, gbw_min, gbw_max)

        # Plot
        comment1 = 'Readout type is ' + dumptype
        comment2 = 'Channel {0:1d}'.format(chan)
        comment3 = 'GainDRE = ' + gain_dre
        comment4 = 'GBW = {0:6.0f} Hz'.format(gbw) 
        fig = plt.figure(figsize=(7, 5))
        min_amp_db, max_amp_db = -80, 10
        ax = fig.add_subplot(1, 1, 1)
        ax.semilogx(f, a_db, '.')
        ax.semilogx([gbw, gbw], [min_amp_db, a_db[i_bw]], '--', color='b', label='Measured bandwidth')
        ax.semilogx([gbw_min, gbw_min], [min_amp_db, max_amp_db], '--', color='r', label='Bandwidth limits')
        ax.semilogx([gbw_max, gbw_max], [min_amp_db, max_amp_db], '--', color='r')
        ax.semilogx([fs_sur2, fs_sur2], [min_amp_db, max_amp_db], '-', color='k', label='Fs / 2')
        ax.legend()
        ax.text(20, min_amp_db + 40, comment1, color='b')
        ax.text(20, min_amp_db + 35, comment2, color='b')
        ax.text(20, min_amp_db + 30, comment3, color='b')
        ax.text(20, min_amp_db + 25, comment4, color=gbw_color)
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

    return(f, a_db)

# -----------------------------------------------------------------------