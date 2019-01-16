# -----------------------------------------------------------------------------
# Imports
import numpy as np
import os
import matplotlib.pyplot as plt
import get_data

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
 
    freqs = file_nb = np.zeros(len(fnames))
    fichlist = fnames.copy()
    for i in range(len(fnames)):
        name=fnames[i]
        file_nb=int(name[i_fich_deb: i_fich_fin])
        fichlist[file_nb]=name 
        freqs[file_nb]=float(name[i_freq_deb:i_freq_fin])
    return(fichlist, freqs)
        
# -----------------------------------------------------------------------
def get_BW(amplitudesdB, leveldB=3):
    r"""
        This function measures the bandwith of a low pass band at a
        specified level.
        
        Parameters
        ----------
        amplitudesdB : array
        contains the values of the amplitudes.

        leveldB : number
        level in dB to define the out of band limit (default = 3dB)

        Returns
        ------- 
        i_max : number
        the index of the last "in-band amplitude"
        """
    
    i_inband = np.where(amplitudesdB >= np.max(amplitudesdB-1*leveldB))
    i_max = i_inband[0][-1]
    return(i_max)

# -----------------------------------------------------------------------
def process_GBW(fulldirname, config, dumptype, chan):
    r"""
        This function reads data from several DRE IQ data files
        to compute the Gain BadWidth product (GBW).
        
        Parameters
        ----------
        dirname : string
        The name of the directory containing the data files
        (no path)

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
    fs = 20e6
    pix_test=40

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    if dumptype == "IQ-ALL":
        pltfilename = os.path.join(plotdirname, "PLOT_GBW_FILTERED.png")
    else:
        pltfilename = os.path.join(plotdirname, "PLOT_GBW_UNFILTERED.png")

    fichlist, f = get_files_freq(datadirname, dumptype)

    nfiles = len(fichlist)
    amplitudes = np.zeros(nfiles)

    print("{0:3d} files to process...".format(nfiles))
    file_nb=0

    for fich in fichlist:
        print('Processing file {0:3d}/{1:3d} '.format(file_nb, nfiles), end='')
        print(fich)
        if dumptype=="IQ-ALL":
            Chan0_i, Chan0_q, Chan1_i, Chan1_q, _ = get_data.readIQ(os.path.join(datadirname, fich))
            if chan==0:
                modulus = \
                    np.sqrt(Chan0_i[:,pix_test].astype('float')**2 + Chan0_q[:,pix_test].astype('float')**2)
            else:
                modulus = \
                    np.sqrt(Chan1_i[:,pix_test].astype('float')**2 + Chan1_q[:,pix_test].astype('float')**2)    
        else: # dumptype="IQ-TEST"
            data, _ = get_data.readfile(os.path.join(datadirname, fich))
            modulus = np.sqrt(data[1:,0].astype('float')**2 + data[1:,1].astype('float')**2)

        amplitudes[file_nb]=np.max(modulus)-np.min(modulus)
                   
        file_nb+=1

    adB = 20*np.log10(amplitudes/np.max(amplitudes))
    # Measuring 3dB BW
    iBW = get_BW(adB, 3)
    GBW = f[iBW]
    GBW_min = 12e3
    GBW_max = 18e3
    FSsur2 = fs/(2*2**7)

    i_gain_deb, i_gain_fin = 46, -4
    gainDRE = fichlist[0][i_gain_deb: i_gain_fin]

    if GBW > GBW_min and GBW < GBW_max:
        GBW_color = 'g'
    else:
        GBW_color = 'r'

    # Plot
    comment1 = 'Readout type is ' + dumptype
    comment2 = 'Channel {0:1d}'.format(chan)
    comment3 = 'GainDRE = ' + gainDRE
    comment4 = 'GBW = {0:6.0f} Hz'.format(GBW) 
    fig = plt.figure(figsize=(7, 5))
    minAmpdB, maxAmpdB = -80, 10
    ax = fig.add_subplot(1, 1, 1)
    ax.semilogx(f, adB, '.')
    ax.semilogx([GBW, GBW], [minAmpdB, adB[iBW]], '--', color='b', label='Measured bandwidth')
    ax.semilogx([GBW_min, GBW_min], [minAmpdB, maxAmpdB], '--', color='r', label='Bandwidth limits')
    ax.semilogx([GBW_max, GBW_max], [minAmpdB, maxAmpdB], '--', color='r')
    ax.semilogx([FSsur2, FSsur2], [minAmpdB, maxAmpdB], '-', color='k', label='Fs / 2')
    ax.legend()
    ax.text(20, minAmpdB + 40, comment1, color='b')
    ax.text(20, minAmpdB + 35, comment2, color='b')
    ax.text(20, minAmpdB + 30, comment3, color='b')
    ax.text(20, minAmpdB + 25, comment4, color=GBW_color)
    ax.axis([10, 1e6, minAmpdB, maxAmpdB])
    ax.set_xlabel(r'Frequency (Hz)')
    ax.set_ylabel(r'Amplitude (dB)')
    ax.set_title(r'Bandwidth measurement')
    major_ticks = np.arange(minAmpdB, maxAmpdB, 10)
    minor_ticks = np.arange(minAmpdB, maxAmpdB, 2)
    ax.set_yticks(major_ticks)
    ax.set_yticks(minor_ticks, minor=True)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)

    fig.tight_layout()
    plt.savefig(pltfilename, bbox_inches='tight')

    return(f, adB)

# -----------------------------------------------------------------------