import numpy as np
from numpy.fft import rfft
import os
import matplotlib.pyplot as plt
import get_data

# -----------------------------------------------------------------------------
def process_dump(fulldirname, config, fs=20e6, Max_duration=0.2):
    r"""
        This function reads and process the data of DRE-DEMUX data dumps.
        
        Parameters
        ----------
        fulldirname : string
        The name of the dump file (with the path)

        config : dictionnary
        Contains path and constants definitions

        fs : number, optional
        Sampling frequency (default is 20 MHz)

        Max_duration : number, optional
        Maximum duration in seconds to be considered for analysis (default is 1)

        Returns
        -------
        Nothing

        """

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    nba, NameA = 12, "INPT" # INPUT signal over 12 bits
    nbb, NameB = 16, "FBCK" # FEEDBACK signal over 16 bits
    nbc, NameC = 16, "BIAS" # BIAS signal over 16 bits

    f_type_deb, f_type_fin = 21, 27
    dumpfilenames1 = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat'\
                and f[f_type_deb:f_type_fin]=="IN-FBK"]
    dumpfilenames2 = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat'\
                and f[f_type_deb:f_type_fin]=="IN-BIA"]

    if len(dumpfilenames1)>0 and len(dumpfilenames2)>0:
        dumpfilename1 = os.path.join(datadirname, dumpfilenames1[0])
        dumpfilename2 = os.path.join(datadirname, dumpfilenames2[0])
        logfilename  = os.path.join(fulldirname, "dumps.log")
        plotfilenameA = os.path.join(plotdirname, "PLOT_DUMP_" + NameA + ".png")
        plotfilenameB = os.path.join(plotdirname, "PLOT_DUMP_" + NameB + ".png")
        plotfilenameC = os.path.join(plotdirname, "PLOT_DUMP_" + NameC + ".png")

        # Getting the data from dump file
        data1, _ = get_data.readfile(dumpfilename1)
        data2, _ = get_data.readfile(dumpfilename2)

        channel=int(data1[0, 1]/2**12)
        a=data1[1:,0]
        b=data1[1:,1]
        c=data2[1:,1]
    
        nval=len(a)
        duration = (nval/fs)
    
        # Reduction of the number of samples (if needed)
        duration = min(duration, Max_duration)
        nval = int(duration * fs)

        # reduction of the number of values to a power of 2 (for a faster FFT)
        power2 = int(np.log(nval)/np.log(2))
        nval = int(2**power2)   

        t = np.linspace(0, duration, nval)

        a = a[:nval]
        b = b[:nval]
        c = c[:nval]

        f_res = fs / nval

        io_str1 = '\n\
    =============================================================================\n\
    =================================  DUMP ANALISYS ============================\n\
    =============================================================================\n\
    Dump file name--------> ' + dumpfilename1 + ' \n\
                            ' + dumpfilename2 + ' \n\
    Log file name---------> ' + logfilename + ' \n\
    Plots saved in files--> ' + plotfilenameA + ' \n\
                            ' + plotfilenameB + ' \n\
                            ' + plotfilenameC + ' \n\
    Channel---------------> {0:2d}\n\
    \n'.format(channel)

        io_str2 = '\
    Number of values in data file----> {0:6d} \n'.format(nval)

        io_str3 = '\
    Duration of measurement----------> {0:6.2f} s \n'.format(duration)

        io_str4 = '\
    =============================================================================\n\
    Dump frequency resolution--------> {0:6.2f} Hz\n\
    Bandwidth factor wrt 1Hz---------> {1:6.2f} db\n\
    =============================================================================\n'\
        .format(f_res, 10*np.log10(f_res))

        flog = open(logfilename, 'w')
        flog.write(io_str1)
        flog.write(io_str2)
        flog.write(io_str3)
        flog.write(io_str4)

        plot_str = ' Number of samples--> {0:8d}   Dump duration--> {1:6.2f} s   Frequency resolution--> {2:8.2f} Hz\n' \
        .format(nval, nval / fs, fs / nval)

    
        ###########################################################################
        print("------------------------------------")
        if abs(a).max()==0:
            print("Data stream A is empty!")
        else:
            print("Processing datastream of INPUT signal...")
            afdb, a_carriers, a_Noise_Power, io_str5 = makeanalysis(a, f_res, nba, config)
            if len(a_carriers) == 0:
                a_carriers = np.ones(1) * 1e6   # if no carriers, a virtual one is 
                                                # defined at 1MHz just to define a 
                                                # focus for the plot
            io_str5 = '\
    == Measurements on data stream ' + NameA + '\n\
    =============================================================================\n'\
            + io_str5
            flog.write(io_str5)
            plot_strA = "Signal "+ NameA + "  " + plot_str
            fmin = (a_carriers[0] - config['ScBandMax']) / 1e6
            fmax = (a_carriers[0] + config['ScBandMax']) / 1e6
            makeplots(t, a, nba, afdb, fs, a_Noise_Power, fmin, fmax, \
                    plot_strA, plotfilenameA)
    
        ###########################################################################
        print("------------------------------------")
        if abs(b).max()==0:
            print("Data stream B is empty!")
        else:
            print("Processing datastream of FEEDBACK signal...")
            bfdb, b_carriers, b_Noise_Power, io_str6 = makeanalysis(b, f_res, nbb, config)
            if len(b_carriers) == 0:
                b_carriers = np.ones(1) * 1e6   # if no carriers, a virtual one is 
                                                # defined at 1MHz just to define a 
                                                # focus for the plot
            io_str6 = '\
    == Measurements on data stream ' + NameB + '\n\
    =============================================================================\n'\
            + io_str6
            flog.write(io_str6)
            plot_strB = "Signal "+ NameB + "  " + plot_str
            fmin = (b_carriers[0] - config['ScBandMax']) / 1e6
            fmax = (b_carriers[0] + config['ScBandMax']) / 1e6
            makeplots(t, b, nbb, bfdb, fs, b_Noise_Power, fmin, fmax, \
                    plot_strB, plotfilenameB)

        ###########################################################################
        print("------------------------------------")
        if abs(c).max()==0:
            print("Data stream C is empty!")
        else:
            print("Processing datastream of BIAS signal...")
            cfdb, c_carriers, c_Noise_Power, io_str7 = makeanalysis(c, f_res, nbc, config)
            if len(c_carriers) == 0:
                c_carriers = np.ones(1) * 1e6   # if no carriers, a virtual one is 
                                                # defined at 1MHz just to define a 
                                                # focus for the plot
            io_str7 = '\
    == Measurements on data stream ' + NameC + '\n\
    =============================================================================\n'\
            + io_str7
            flog.write(io_str7)
            plot_strC = "Signal "+ NameC + "  " + plot_str
            fmin = (c_carriers[0] - config['ScBandMax']) / 1e6
            fmax = (c_carriers[0] + config['ScBandMax']) / 1e6
            makeplots(t, c, nbc, cfdb, fs, c_Noise_Power, fmin, fmax, \
                    plot_strC, plotfilenameC)
    
        flog.close()

    print("Done!")
    print("------------------------------------")
# -----------------------------------------------------------------------------

def makeanalysis(sig, df, nb, config):
    r"""
        This function does the analysis of a signal. In time domain (measurment
        of the power) but mostly in frequency domain (fft, normalisation, ...).

        Parameters
        ----------
        sig : array_like
        The (time domain) signal to be analysed.

        df : number
        The spectral resolution of the signal.

        nb : number
        Number of bits used to code the digital signal (16 for the DAC).

        config : dictionnary
        Contains path and constants definitions

        Returns
        -------
        sigfdb : array_like
        The signal in frequency domain and expressed in dB FSR.

        carriers : array_like
        Indexes of detected peaks according to the input array

        Noise_Power : number
        The maximum power measured in the science band around the carriers.

        io_str : string
        The log message.

    """
    
    nval=len(sig)
    #######################################################################
    # Computation of fft
    sigf = abs(rfft(sig))
     
    #######################################################################
    # Get the peak-to-peak amplitude of the time signal.
    PeakPeak = 2.*max(abs(sig))

    # FSR power level (constant signal equal to DAC at full scale)
    Efsr_db = 10*np.log10(nval*(2**nb)**2)

    # Power of the signal measured in time domain
    Et_db = 10*np.log10(np.sum(sig.astype(float)**2))
    # Power obtained in frequency domain (to be corrected for)
    Ef_db = 10*np.log10(np.sum(sigf.astype(float)**2))

    sigfdb = -300.*np.ones(len(sigf))
    inotzero = np.where(sigf != 0)[0]
    sigfdb[inotzero] = 20*np.log10(sigf[inotzero]) - Ef_db + Et_db - Efsr_db

    ##### Noise power measurement around each carriers
    fcarriers=np.array([])
    Noise_Power=0
    if nb > 12: # bias or feedback
        print("Peak detection...")
        icarriers = peakdetect(sigfdb)
        fcarriers = icarriers*df
        ncar=len(fcarriers)
        print(" => {0:3d} carriers detected. \nCarrier frequencies (Hz):".format(ncar))
        freq_str, power_str =  " > ",  " > "
        if ncar > 0 and ncar < 80:
            for carrier in range(ncar):
                freq_str += '{0:12.2f}'.format(fcarriers[carrier])
                if (carrier+1) % 8 == 0:
                    freq_str += "\n > "
            print(freq_str)
            print("\nNoise power calculation... \nMeasured power (dBFS / 2 x Binfo)")
            Noise_Power = noiseandspurpower(sigf, icarriers, df, config) - Ef_db + Et_db - Efsr_db
            for carrier in range(ncar):
                power_str += '{0:12.2f}'.format(Noise_Power[carrier])
                if (carrier+1) % 8 == 0:
                    power_str += "\n > "
            print(power_str)
    else:
        print("Calculation of mean noise power... (dBFS / 2 x Binfo)")
        Noise_Power = 10*np.log10(np.sum(sig**2)) - 10*np.log10(df*nval/2) - Efsr_db 
        print(Noise_Power)
 
    io_str = '\
Peak Peak amplitude--------------> {0:9.2f} \n\
Measured crest factor------------> {1:9.2f} \n\
Max of spectrum------------------> {2:9.2f} dB\n'\
    .format(PeakPeak, crestfactor(sig), sigfdb.max())
    if nb > 12:
        io_str += 'Number of carriers detected------> {0:9d} : \n'\
        .format(ncar)
        io_str += freq_str + '\n\
Noise power around carriers-> (dBFS / 2 x Binfo)\n\
=============================================================================\n'\
        + power_str

    return(sigfdb, fcarriers, Noise_Power, io_str)


# -----------------------------------------------------------------------------
def makeplots(t, sig, nb, sigfdb, fs, Noise_Power, fmin, fmax, io_str, plotfilename):
    r"""
        This function plots a signal in time and frequency domain.

        Parameters
        ----------
        t : array_like
        The time series (x values of time domain plots)

        sig : array_like
        The data in time domain.

        nb : number
        The number of bits used to code the signal.

        sigfdb : array_like
        The data in frequency domain (expressed in dbFSR)

        fs : number
        The sampling frequency of the data.

        Noise_Power : number
        The maximum power measured in the science band around the carriers.

        fmin / fmax : numbers
        The limits of the x range for the frequency domain plots.

        io_str : string
        text to be displayed on the plot

        plotfilename : string
        name of the plotfile.

        Returns
        -------
        Nothing

    """
    Cf = crestfactor(sig-sig.mean())
    moyenne = sig.mean()
    moyenne_dBFS = 20*np.log10(2**nb/np.abs(moyenne))
    PeakPeak = 2.*abs(sig).max()
    FSR_over_PeakPeak = 2**nb/PeakPeak
    io_str2 = '\n Crest factor = {0:5.2f}    FSR_DAC/PeakPeak = {1:6.2f}\n Mean = {2:5.2f} = {3:5.2f} dBFS'\
        .format(Cf, FSR_over_PeakPeak, moyenne, moyenne_dBFS)

    f = np.linspace(0, fs/2, len(sigfdb))
    L1 = 1024         # length of time plot (long)
    fig = plt.figure(figsize=(11, 8))
    fig.text(0.05, 0.98, io_str, family='monospace')
    ytimetitle = "amplitude"
    yfreqtitle = "amplitude (dB FSR)"

    # some reference levels
    AllSinesRmsLevel = -20.*np.log10(2*np.sqrt(2.)*4.5*np.sqrt(40))
    # some reference levels
    NoiseReq = -156

    # setting of the plot-ranges
    fxMin2 = 0.7
    fxMax2 = 5.2
    fyMin = -180
    fyMax = 0

    #####################################################################
    # Plotting data in time domain

    ax1 = fig.add_subplot(2, 2, 1)
    #ax1.plot(1000*t, sig, 'b')
    #ax1.plot(1000*t, sig, '.r')
    ax1.plot(1000*t[0:L1], sig[0:L1], 'b')
    ax1.text(0, 2**(nb-1)*0.75, io_str2, color='b')
    #ax1.plot(1000*t[0:L1], sig[0:L1], '.r')
    ax1.set_ylabel(ytimetitle)
    ax1.set_xlabel('Time (ms)')
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    #ax1.axis([0, 1000*t[-1], -1*2**(nb-1), 2**(nb-1)])
    ax1.axis([0, 1000*L1/fs, -1*2**(nb-1), 2**(nb-1)])

    #####################################################################
    # Plotting data in frequency domain

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(f/1e6, sigfdb, 'b', linewidth=1)
    if nb==16:
        ax2.plot([0, f[-1]], [AllSinesRmsLevel, AllSinesRmsLevel], '--g', linewidth=1.5)
        ax2.plot([0, f[-1]], [NoiseReq, NoiseReq], '--k', linewidth=1.5)
        #ax2.plot([0, f[-1]], [Noise_Power, Noise_Power], '--k', linewidth=1.5)
    ax2.set_ylabel(yfreqtitle)
    ax2.set_xlabel('Frequency (MHz)')
    ax2.grid(color='k', linestyle=':', linewidth=0.5)
    ax2.axis([0, fs/2 / 1e6, fyMin, fyMax])

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(f / 1e6, sigfdb, 'b', linewidth=1)
    #ax3.plot(f / 1e6, sigfdb, '.r')
    if nb==16:
        ax3.plot([0, f[-1]], [AllSinesRmsLevel, AllSinesRmsLevel], '--g', linewidth=1.5)
        ax3.plot([0, f[-1]], [NoiseReq, NoiseReq], '--k', linewidth=1.5)
    #ax3.plot([0, f[-1]], [Noise_Power, Noise_Power], '--k', linewidth=1.5)
    ax3.set_ylabel(yfreqtitle)
    ax3.set_xlabel('Frequency (MHz)')
    ax3.grid(color='k', linestyle=':', linewidth=0.5)
    ax3.axis([fxMin2, fxMax2, fyMin, fyMax])

    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(f / 1e6, sigfdb, 'b', linewidth=1)
    ax4.plot(f / 1e6, sigfdb, '.r')
    if nb==16:
        ax4.plot([0, f[-1]], [AllSinesRmsLevel, AllSinesRmsLevel], '--g', linewidth=1.5)
        ax4.plot([0, f[-1]], [NoiseReq, NoiseReq], '--k', linewidth=1.5)
    #ax4.plot([0, f[-1]], [Noise_Power, Noise_Power], '--k', linewidth=1.5)
    ax4.set_ylabel(yfreqtitle)
    ax4.set_xlabel('Frequency (MHz)')
    ax4.grid(color='k', linestyle=':', linewidth=0.5)
    ax4.axis([fmin, fmax, fyMin, fyMax])

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename, bbox_inches='tight')
    print('Plots done')

# -----------------------------------------------------------------------------
def peakdetect(sig, margin=6):
    r"""
        This function detects all peaks in a 1D signal. A threshold is used
        to define the lower limit (wrt the maximum of the signal) for the search.

        Parameters
        ----------
        sig : array_like
        The signal in which peaks will be detected

        margin : number
        The range where to look for maxima wrt to the absolute maximum 
        (default is 6).

        Returns
        -------
        detected_peaks_indexes : array_like
        Indexes of detected peaks according to the input array

        Example
        -------
        >>> a = peakdetect(vector,margin)
        >>> a
        array([ 5,  9,  23, ...,  102,  143,  340])

        """
    # First look for maximum areas
    ithreshold = np.where(sig > sig.max() - margin)[0]

    # Then look for local maxima in each area
    derivative1 = sig[1:] - sig[:-1]
    i_zero = np.where(derivative1 == 0)[0]       # To avoid division
    derivative1[i_zero] = 1                      # by 0 at the next step.
    derivative2 = derivative1 / abs(derivative1) # keeps only +1 and -1
    pic = derivative2[:-1] - derivative2[1:]     # equals 2 at local maxima
    i = np.where(pic[ithreshold-1] == 2)[0]
    return(ithreshold[i])

# -----------------------------------------------------------------------------
def crestfactor(signal):
    r"""
        This function computes the crest factor (cf) of a signal.

        The cf is defined as the ratio between the peak value of a signal and
        its rms value.

        Parameters
        ----------
        signal : array_like
        an array

        Returns
        -------
        cf : number
        The crest factor computed using cf = Peak / rms

        Example
        -------
        >>> crestfactor(array([3., 7., 1., 5., 4., 6., 3.5]))
        0.651887905613

        """
    peak = abs(signal).max()
    return(peak / signal.std())

# -----------------------------------------------------------------------------
def noiseandspurpower(sigf, indexes, df, config):
    r"""
        This function measures the power of the signal in a band around
        specific frequency indexes.

        Parameters
        ----------
        sig : array_like
        The signal in which the power is measured.

        indexes : array_like
        An array containing the indexes of the carriers.

        df : number
        The frequency resolution.

        config : dictionnary
        Contains path and constants definitions

        Returns
        -------
        The values of the power measured around each carrier.

    """

    Powers = np.zeros(len(indexes))
    for k, index in enumerate(indexes):
        l_side_1 = index - int(np.ceil(config['ScBandMax']/df))
        l_side_2 = index - int(np.ceil(config['ScBandMin']/df))
        r_side_1 = index + int(np.ceil(config['ScBandMin']/df))
        r_side_2 = index + int(np.ceil(config['ScBandMax']/df))
        #npts = l_side_2-l_side_1 + r_side_2-r_side_1
        Powers[k] = np.sum(sigf.astype(float)[l_side_1:l_side_2]**2) \
                    + np.sum(sigf.astype(float)[r_side_1:r_side_2]**2)

    PowersdB=-1*np.inf*np.ones(len(Powers))
    inotzero=np.where(Powers != 0)[0]
    PowersdB[inotzero] = 10*np.log10(Powers[inotzero])
    return(PowersdB)

# -----------------------------------------------------------------------
def process_dump_pulses(fulldirname, config, fs=20e6, Max_duration=0.2):
    r"""
        This function reads and process the data of DRE-DEMUX data dumps.
        INPUT signal while pulses.
        
        Parameters
        ----------
        fulldirname : string
        The name of the dump file (with the path)

        config : dictionnary
        Contains path and constants definitions

        fs : number, optional
        Sampling frequency (default is 20 MHz)

        Max_duration : number, optional
        Maximum duration in seconds to be considered for analysis (default is 1)

        Returns
        -------
        Nothing

        """

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    nb = 12 # INPUT signal over 12 bits

    f_type_deb, f_type_fin = 21, 33
    dumpfilenames = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat'\
                and f[f_type_deb:f_type_fin]=="IN-BIA_PULSE"]

    if len(dumpfilenames)>0:
        dumpfilename = os.path.join(datadirname, dumpfilenames[0])
        plotfilename = os.path.join(plotdirname, "PLOT_INPUT_WITH_PULSES.png")

        # Getting the data from dump file
        data, _ = get_data.readfile(dumpfilename)

        channel=int(data[0, 1]/2**12)
        a=data[1:,0]
        pp=a.max()-a.min()
    
        nval=len(a)
        duration = (nval/fs)
    
        t = np.linspace(0, duration, nval)

        Pulse_length = 1024*128
        center = a[Pulse_length:-1*Pulse_length] # Only consider the center of the dumb to have a full pulse
        iMax = np.where(center==np.max(center))[0][0]+Pulse_length
        ideb = iMax - 128*100
        ifin = iMax + Pulse_length

        # Making plot
        fig = plt.figure(figsize=(11, 8))
        io_str="File: "+dumpfilenames[0]+", Channel {0:1d}, FSR(ADC)/PeakPeak = {1:5.2f}".format(channel, 2.**12/pp)
        fig.text(0.1, 0.982, io_str, family='monospace')

        ax = fig.add_subplot(1, 1, 1)
        ax.plot(t[ideb:ifin]*1e3, a[ideb:ifin])
        ax.set_ylabel("ADC unit (FSR range)")
        ax.set_xlabel("Time (ms)")
        ax.grid(color='k', linestyle=':', linewidth=0.5)
        ax.set_ylim([-2**(nb-1), 2**(nb-1)])

        fig.tight_layout()
        #plt.show()
        plt.savefig(plotfilename, bbox_inches='tight')

# -----------------------------------------------------------------------
def plot_spectra(sptdB, fs, config, pltfilename, Cf, FSR_over_PeakPeak, ncar, suffixe, BW_correction_factor_dB, pix_zoom=40):
    r"""
        This function plots spectra of IQ data wrt the carriers.
                
        Parameters
        ----------
        sptdB : arrays
        Spectra.
        
        fs: number
        sampling rate of the data (default = 20e6/2**7 = 156.25 kHz)

        config: dictionnary
        contains usefull parameters
        
        pltfilename, suffixe: strings
        Used to specify the plot file names.

        Cf : number.
        Crest factor of the signal (used for the conversion dBc / dBFS)

        FSR_over_PeakPeak : number.
        Ratio between signal peak to peak value and DAC FSR (used for the conversion dBc / dBFS)

        ncar : number.
        Number of carriers

        suffixe : string. Extension to be added to the plot file name.

        BW_correction_factor_dB : number.
        BW correction factor that has been applied to obtain the spectra in dB/Hz. This value needs to be 
        indicated on the plot so that it can be taken into account when measuring spuriouses.

        pix_zoom: number
        The id of a pixel to be zoomed in.
        (default value is 40 = test pixel)

        Returns
        -------
        Nothing           
        """

    npts = len(sptdB[0,:])
    f=np.arange(npts)*fs/(2*npts)
    Impact_2_DACs = 3    # Because the test set-up involves two DACs (bias and feedback)
    Impact_non_stationarity = 20*np.log10(1.4) # Because our test set-up does not reproduce the TES non-stationarity 
    Impact_BBFB = 3  # BBFB adds 3 dB noise (TBC)
    SNR_pix_min = config['SNR_pix'] + Impact_2_DACs + Impact_non_stationarity
    SNR_pix_max = SNR_pix_min + Impact_BBFB
    ScBandMin = config['ScBandMin'] 
    ScBandMax = config['ScBandMax'] 
    BW_res_warn = r'A BW res. correction factor of {0:4.2f}dB has been applied on the spectra. It shall be corrected for spurious measurements.'.format(BW_correction_factor_dB)

    # Plot for zoom pixel only
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.semilogx(f[1:], sptdB[pix_zoom,1:])
    ax.semilogx([ScBandMin, ScBandMax], [SNR_pix_min, SNR_pix_min], ':r', linewidth=3)
    ax.semilogx([ScBandMin, ScBandMax], [SNR_pix_max, SNR_pix_max], ':r', linewidth=3)
    ax.semilogx([ScBandMin, ScBandMin], [SNR_pix_min, 0], ':r', linewidth=3)
    ax.semilogx([ScBandMax, ScBandMax], [SNR_pix_min, 0], ':r', linewidth=3)
    ax.text(40, -17, r'band of interest', fontsize=11)
    ax.text(ScBandMin, SNR_pix_min-5, r'DRD requirement level', fontsize=11)
    ax.text(1.5, -159, BW_res_warn, fontsize=11)
    L1, L2=100, 2.5
    ax.arrow(ScBandMin, -20, ScBandMax-L1, 0, head_width=3, head_length=L1, fc='k', ec='k')
    ax.arrow(ScBandMax, -20, -ScBandMax+ScBandMin+L2, 0, head_width=3, head_length=L2, fc='k', ec='k')
    ax.axis([1, f[-1], -160, 0])
    ax.set_xlabel(r'Frequency (Hz)')
    ax.set_ylabel(r'DRD (dBc/Hz)')
    ax.set_title(r'Pixel {0:2d}'.format(pix_zoom))
    # Major ticks every 20, minor ticks every 10
    major_ticks = np.arange(-160, 1, 20)
    minor_ticks = np.arange(-160, 1, 10)
    #ax.set_xticks(major_ticks)
    #ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks)
    ax.set_yticks(minor_ticks, minor=True)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_weight('bold')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(17)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    if Cf!=0:
        Carrier_vs_FS_dB=20*np.log10(FSR_over_PeakPeak*2*Cf*np.sqrt(ncar))
        ax2 = plt.gca().twinx()
        ax2.axis([f[1], f[-1], -160-Carrier_vs_FS_dB, 0-Carrier_vs_FS_dB])
        ax2.set_ylabel(r'DRD (dBFS/Hz)')
        ax2.yaxis.label.set_weight('bold')
        ax2.yaxis.label.set_fontsize(17)
        for item in (ax2.get_yticklabels()):
            item.set_fontsize(15)

    fig.tight_layout()
    plt.savefig(pltfilename+suffixe+'_zoom'+str(pix_zoom)+'.png', bbox_inches='tight')

    # Plot for all other pixels (but the test pixels)
    n_boxes=40 
    n_lines=5
    n_cols =8

    # Checking which pixel is on
    pix_ON = np.ones((n_boxes), dtype=bool)
    for pix in range(n_boxes):
        if sptdB[pix,1:].max()==sptdB[pix,1:].min():
            pix_ON[pix]=False
            print("\n------------------> Pixel {0:2d} is off.".format(pix))

    fig = plt.figure(figsize=(18, 12))
    for box in range(n_boxes):
        if pix_ON[box]:
            ax = fig.add_subplot(n_lines, n_cols, box+1)
            ax.semilogx(f[1:], sptdB[box,1:])
            ax.semilogx([ScBandMin, ScBandMax], [SNR_pix_min, SNR_pix_min], ':r')
            ax.semilogx([ScBandMin, ScBandMax], [SNR_pix_max, SNR_pix_max], ':r')
            ax.semilogx([ScBandMin, ScBandMin], [SNR_pix_min, 0], ':r')
            ax.semilogx([ScBandMax, ScBandMax], [SNR_pix_min, 0], ':r')
            ax.axis([f[1], f[-1], -160, 0])
            ax.grid(color='k', linestyle=':', linewidth=0.5)
            ax.set_title(r'Pixel {0:2d}'.format(box))
            if box//n_cols == n_lines-1:
                ax.set_xlabel(r'Frequency (Hz)')
            else:
                plt.xticks(visible=False)
            if box%n_cols == 0:
                ax.set_ylabel(r'DRD (dBc/Hz)')
            else:
                plt.yticks(visible=False)
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
                item.set_fontsize(15)
            for item in (ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(13)
    fig.tight_layout()
    plt.savefig(pltfilename+suffixe+'_40pix.png', bbox_inches='tight')

# -----------------------------------------------------------------------
def get_Cf_and_FSRoverPeakPeak_from_file(fulldirname, quiet=True):
    r"""   
        This function gets data from a ADC - Feedback dump file, it 
        computes the crest factor and the FSR_over_peakpeak ratio.

        Parameters
        ----------
        fulldirname : string
        Name of the directory
        """

    f_type_deb, f_type_fin = 21, 27
    dumpfilenames = [f for f in os.listdir(fulldirname) \
            if os.path.isfile(os.path.join(fulldirname, f)) \
            and f[-4:]=='.dat'\
            and f[f_type_deb:f_type_fin]=="IN-BIA"]

    if len(dumpfilenames)>0:
        filename = os.path.join(fulldirname, dumpfilenames[0])

        data, _ = get_data.readfile(filename)
        feedback = data[1:,1].astype(float)
        Cf = crestfactor(feedback)
        PeakPeak = 2.*max(abs(feedback))
        FSR_over_PeakPeak = 2.**16/PeakPeak
        
        # Searching the carriers
        feedback_F = abs(rfft(feedback[:2**16]))
        feedback_Fdb = -300.*np.ones(2**16)
        inotzero = np.where(feedback_F != 0)[0]
        feedback_Fdb[inotzero] = 20*np.log10(feedback_F[inotzero])
        ncar = len(peakdetect(feedback_Fdb))
        #print("============>>>>>>>>>>>> ", ncar)

        if not quiet:
            print("Feedback crest factor: {0:5.2f}".format(Cf))
            print("Feedback peak value:   {0:5.2f}".format(PeakPeak))
            print("Number of carriers:    {0:5d}".format(ncar))
            print("DAC peak value:        {0:5d}".format(2**16-1))
            print("FSR / Peak:            {0:5.2f}".format(FSR_over_PeakPeak))
    else:
        Cf = FSR_over_PeakPeak = ncar = 0

    return(Cf, FSR_over_PeakPeak, ncar)
 
# -----------------------------------------------------------------------
def processIQ_multi(fulldirname, config, fs=20e6, pix_zoom=40, window=False, BW_CORRECTION=True):
    r"""
        This function reads data from several DRE IQ data files.
        It computes the accumulated spectra and makes the plots.
        
        Parameters
        ----------
        dirname : string
        The name of the directory containing the data files
        (no path)
        
        fs: number
        The data sampling rate at DAC and ADC .
        (default value is 20e6)
        
        pix_zoom: number
        The id of a pixel to be zoomed in.
        (default value is 40 = test pixel)

        window: boolean
        Specifies if a windowing shall be done before the FFT.
        If True a Blackman Harris window is applied. (Default is False)

        BW_CORRECTION: boolean
        Specifies if a correction factor has to be applied to take into account the resolution BW (Default is True).
        If this factor is applied the measurment will be wrong for spuriouses (Default is True).

        Returns
        ------- 
        spt0dB : Array containing the accumulated pixels spectra
           
        """

    # computing sampling frequency at BBFB output
    fs = fs/2**7

    npix=41
    npts_max=2**17
    nb_short_files=0

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_carrier_spt")

    i_test_deb, i_test_fin = 21, 40
    test = "IQ-ALL_Science-Data"
    fichlist = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat' \
                and f[i_test_deb:i_test_fin]==test]

    if len(fichlist)>0:
        print('Measurment of signal crest factor from feedback dump files if they exist')
        Cf0, FSR_over_PeakPeak0, ncar = get_Cf_and_FSRoverPeakPeak_from_file(datadirname)

        # definning file length
        Chan0_i, Chan0_q, _, _, FLAG_ERROR = get_data.readIQ(os.path.join(datadirname, fichlist[0]))
        npts = len(Chan0_i[:,0])
        print('Npts:', npts)
        npts=min(npts_max, 2**np.int(np.log(npts)/np.log(2)))
        duration=npts/fs
        # Factor to correct the resolution BW effect on noise
        BW_correction_factor_dB=10*np.log10(duration)
        print("Scan duration is about: {0:6.4f}".format(duration))
        print("WARNING, a bandwidth correction factor of {0:6.4f}dB is applied on the spectra.".format(BW_correction_factor_dB))
        print("This offset needs to be taken into account when considering spurious values.")

        spt0dB=-300*np.ones((npix, npts//2+1))
        total_spt0=np.zeros((npix, npts//2+1))
        spt = np.zeros((npix, npts//2+1))

        win = 1
        if window:
            win = np.blackman(npts)
    
        nfiles = len(fichlist)
        print("{0:3d} files to process...".format(nfiles))
        file=0
        errors_counter = 0
        CHAN0_EMPTY=True
        for fich in fichlist:
            file+=1
            print('Processing file {0:3d}/{1:3d} '.format(file, nfiles), end='')
            print(fich)
            Chan0_i, Chan0_q, _, _, FLAG_ERROR = get_data.readIQ(os.path.join(datadirname, fich))
            npts_current = len(Chan0_i[:,0])
            if FLAG_ERROR:
                errors_counter += 1

            if npts_current >= npts:

                Chan0_modulus = \
                    np.sqrt(Chan0_i[0:npts,:].astype('float')**2 + Chan0_q[0:npts,:].astype('float')**2)
                            
                if Chan0_modulus.max() > 0: # data exists
                    CHAN0_EMPTY=False
                    for pix in range(npix):
                        spt[pix,:] = abs(rfft((Chan0_modulus[:,pix])*win))
                    total_spt0 += spt**2
                                
            else:
                print("File is too short!")
                nb_short_files += 1 

        print("Data processing is done.")
        print("{0:4d} corrupted files found.".format(errors_counter))
        print("{0:4d} files were too short for processing.".format(nb_short_files))
        print("Doing the plots...", end='')
    
        if not CHAN0_EMPTY:
            # Normalisation wrt each carrier
            for pix in range(npix): 
                i_good = np.where(total_spt0[pix,:] > 0)[0]
                spt0dB[pix,i_good] = 10*np.log10(total_spt0[pix,i_good])
                spt0dB[pix,:] = spt0dB[pix,:] - spt0dB[pix,:].max()
            # 3 dB correction to compensate the impact of the rfft on DC bin
            spt0dB[:,1:] += 3
            # Normalisation to a RBW of 1Hz
            if BW_CORRECTION:
                spt0dB[:,1:] += BW_correction_factor_dB
            plot_spectra(spt0dB, fs, config, pltfilename, Cf0, FSR_over_PeakPeak0, ncar, '', BW_correction_factor_dB, pix_zoom)
    else:
        spt0dB=0
        
    return(spt0dB)

# -----------------------------------------------------------------------
def processIQ_TST_multi(fulldirname, config, fs=20e6, window=False, BW_CORRECTION=True):
    r"""
        This function reads data from several DRE IQ-TST data files.
        It computes the accumulated spectra and makes the plots.
        
        Parameters
        ----------
        dirname : string
        The name of the directory containing the data files
        (no path)
        
        fs: number
        The data sampling rate at DAC and ADC .
        (default value is 20e6)
        
        window: boolean
        Specifies if a windowing shall be done before the FFT.
        If True a Blackman Harris window is applied. (Default is False)

        BW_CORRECTION: boolean
        Specifies if a correction factor has to be applied to take into account the resolution BW (Default is True).
        If this factor is applied the measurment will be wrong for spuriouses (Default is True).
        
        Returns
        ------- 
        spt0dB : Array containing the accumulated pixels spectra
           
        """

    npts_max=2**26
    nb_short_files=0

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_carrier-TST_spt")

    i_test_deb, i_test_fin = 21, 40
    test = "IQ-TST_Science-Data"
    fichlist = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat' \
                and f[i_test_deb:i_test_fin]==test]

    if len(fichlist)>0:
        print('Measurment of signal crest factor from feedback dump files if they exist')
        Cf, FSR_over_PeakPeak, ncar = get_Cf_and_FSRoverPeakPeak_from_file(datadirname)

        # definning file length
        data, dumptype = get_data.readfile(os.path.join(datadirname, fichlist[0]))
        if dumptype != 9:
            raise ValueError('Wrong dumptype')

        # decommutation des donnees        
        npts = len(data[1:,0])   
        print('Npts:', npts)
        npts=min(npts_max, 2**np.int(np.log(npts)/np.log(2)))
        duration=npts/fs
        # Factor to correct the resolution BW effect on noise
        BW_correction_factor_dB=10*np.log10(duration)
        print("Scan duration is about: {0:6.4f}".format(duration))
        print("WARNING, a bandwidth correction factor of {0:6.4f}dB is applied on the spectra.".format(BW_correction_factor_dB))
        print("This offset needs to be taken into account when considering spurious values.")

        sptdB=np.zeros((npts//2+1))
        total_spt=np.zeros((npts//2+1))

        win = 1
        if window:
            win = np.blackman(npts)
    
        nfiles = len(fichlist)
        print("{0:3d} files to process...".format(nfiles))
        file=0
        errors_counter = 0
        EMPTY=True
        for fich in fichlist:
            file+=1
            print('Processing file {0:3d}/{1:3d} '.format(file, nfiles), end='')
            print(fich)

            data, dumptype = get_data.readfile(os.path.join(datadirname, fich))
            if dumptype != 9:
                raise ValueError('Wrong dumptype')

            # decommutation des donnees        
            npts_current = len(data[1:,0])

            if npts_current >= npts:

                modulus = \
                    np.sqrt(data[1:npts+1,0].astype('float')**2 + data[1:npts+1,1].astype('float')**2)
                            
                if modulus.max() > 0: # data exists
                    EMPTY=False
                    total_spt += abs(rfft(modulus*win))**2
                                
            else:
                print("File is too short!")
                nb_short_files += 1 

        print("Data processing is done.")
        print("{0:4d} corrupted files found.".format(errors_counter))
        print("{0:4d} files were too short for processing.".format(nb_short_files))
        print("Doing the plots...", end='')
    
        if not EMPTY:            
            sptdB = 10*np.log10(total_spt)
            # Normalisation wrt carrier
            sptdB -= sptdB.max()
            # 3 dB correction to compensate the impact of the rfft on DC bin
            sptdB[1:] += 3
            # Normalisation to a RBW of 1Hz
            if BW_CORRECTION:
                sptdB[1:] += BW_correction_factor_dB

            npts = len(sptdB)
            f=np.arange(npts)*fs/(2*npts)
            Impact_2_DACs = 3    # Because the test set-up involves two DACs (bias and feedback)
            Impact_non_stationarity = 20*np.log10(1.4) # Because our test set-up does not reproduce the TES non-stationarity 
            Impact_BBFB = 3  # BBFB adds 3 dB noise (TBC)
            SNR_pix_min = config['SNR_pix'] + Impact_2_DACs + Impact_non_stationarity
            SNR_pix_max = SNR_pix_min + Impact_BBFB
            ScBandMin = config['ScBandMin'] 
            ScBandMax = config['ScBandMax'] 
            BW_res_warn = r'A BW res. correction factor of {0:4.2f}dB has been applied on the spectra. It shall be corrected for spurious measurements.'.format(BW_correction_factor_dB)

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(1, 1, 1)
            ax.semilogx(f[1:], sptdB[1:])
            ax.semilogx([ScBandMin, ScBandMax], [SNR_pix_min, SNR_pix_min], ':r', linewidth=3)
            ax.semilogx([ScBandMin, ScBandMax], [SNR_pix_max, SNR_pix_max], ':r', linewidth=3)
            ax.semilogx([ScBandMin, ScBandMin], [SNR_pix_min, 0], ':r', linewidth=3)
            ax.semilogx([ScBandMax, ScBandMax], [SNR_pix_min, 0], ':r', linewidth=3)
            ax.text(40, -17, r'band of interest', fontsize=11)
            ax.text(ScBandMin, SNR_pix_min-5, r'DRD requirement level', fontsize=11)
            ax.text(1.5, -159, BW_res_warn, fontsize=11)
            L1, L2=100, 2.5
            ax.arrow(ScBandMin, -20, ScBandMax-L1, 0, head_width=3, head_length=L1, fc='k', ec='k')
            ax.arrow(ScBandMax, -20, -ScBandMax+ScBandMin+L2, 0, head_width=3, head_length=L2, fc='k', ec='k')
            #ax.axis([1, f[-1], -180, 0])
            ax.set_ylim(-180, 0)
            ax.set_xlim(1, f[-1])
            ax.set_xlabel(r'Frequency (Hz)')
            ax.set_ylabel(r'DRD (dBc/Hz)')
            ax.set_title(r'Test pixel')
            # Major ticks every 20, minor ticks every 10
            major_ticks = np.arange(-160, 1, 20)
            minor_ticks = np.arange(-160, 1, 10)
            #ax.set_xticks(major_ticks)
            #ax.set_xticks(minor_ticks, minor=True)
            ax.set_yticks(major_ticks)
            ax.set_yticks(minor_ticks, minor=True)
            ax.grid(which='minor', alpha=0.2)
            ax.grid(which='major', alpha=0.5)
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
                item.set_weight('bold')
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
                item.set_fontsize(17)
            for item in (ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(15)
            if Cf!=0:
                Carrier_vs_FS_dB=20*np.log10(FSR_over_PeakPeak*2*Cf*np.sqrt(ncar))
                ax2 = plt.gca().twinx()
                ax2.axis([f[1], f[-1], -160-Carrier_vs_FS_dB, 0-Carrier_vs_FS_dB])
                ax2.set_ylabel(r'DRD (dBFS/Hz)')
                ax2.yaxis.label.set_weight('bold')
                ax2.yaxis.label.set_fontsize(17)
                for item in (ax2.get_yticklabels()):
                    item.set_fontsize(15)

            fig.tight_layout()
            plt.savefig(pltfilename+'.png', bbox_inches='tight')

    else:
        sptdB=0
        
    return(sptdB)

# -----------------------------------------------------------------------
