'''
Created on 26 juil. 2019

@author: philippepeille
'''
import os
import numpy as np
import argparse
import general_tools
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize.minpack import curve_fit

DRE_DATA_ROOT = os.path.normcase("./EP_test_data/") # Put path to folder containing data files here
PIXEL_DATA_DIR = os.path.normcase("./Pixel_data_LPA75um_AR0.5/")
TES_NOISE_SPECTRUM = os.path.join(PIXEL_DATA_DIR,"noise_spectra_bbfb_noFBDAC.fits")
XIFUSIM_TEMPLATE = os.path.join(PIXEL_DATA_DIR,"pulse_withBBFB.npy")
THRESHOLD_ESTIMATE=8000
PREBUFF=200
JITTER_MARGIN=10

# ############################################################
# Function to read data records
# ############################################################
def get_records(filename, verbose=False):
    '''Reads sample records from a file

    Arguments:
        - filename: name of the file containing the sample records
        - verbose: option to print the nonlinear factor

    Returns: an array containing the records
    '''

    # getting record length from file
    dtini=np.dtype([('timestamp', np.float), \
                 ('channelId', np.int8), \
                 ('pixelId', np.int8), \
                 ('recordSize', np.int16)])

    fdat=open(filename, 'rb')
    data=np.fromfile(fdat, dtype=dtini)
    fdat.close()

    recordsize=data[0]['recordSize']    

    # getting the data
    dt=np.dtype([('timestamp', np.float), \
                 ('channelId', np.int8), \
                 ('pixelId', np.int8), \
                 ('recordSize', np.int16), \
                 ('records', np.dtype((np.float, recordsize)))])
    
    fdat=open(filename, 'rb')
    data=np.fromfile(fdat, dtype=dt)
    fdat.close()

    print(data['pixelId'][0])
    print(data['pixelId'][1])
    print(data['records'][0][0])
    print(data['records'][0][1])
    print(data['records'][0][2])
    plt.plot(data['timestamp'])
    #plt.plot(data['records'][0])
    plt.show()

    return data['records']

# ############################################################
# Function to get the pixel non linearity factor from its information file
# ############################################################
def get_nonlinear_factor(dirname, verbose=False):
    '''Reads the pixel's nonlinear factor from its information file

    Arguments:
        - dirname: name of the directory containing the pixel's related data
        - verbose: option to print the nonlinear factor

    Returns: the nonlinear factor
    '''
    filename=os.path.join(dirname,'pixel_info.csv')
    NLF=general_tools.get_csv(filename)['NONLINEAR_FACTOR']
    if verbose:
        print("Pixel non linearity factor: ", NLF)
    return NLF

# ############################################################
# Function to create pulse average
# ############################################################
def pulse_average(pulse_list,ignore_outlayers=False):
    """Creates pulse template from a set of pulses. Outlayers can be manually rejected.
    
    Argumnents:
        - pulse_list: 
        
    Returns: pulse_template 
        - pulse_template: template created from the average of detected pulses
    """
    satisfied=False
    
    # Iteratively compute the average of the detected pulses and reject 10% worst if requested
    pulse_list = np.array(pulse_list)
    while not satisfied:
        mean_pulse = np.mean(pulse_list,0)
        if ignore_outlayers:
            break
        diff_list = np.sum(abs(pulse_list-mean_pulse),1)
        plt.clf()
        plt.plot(mean_pulse,label='Averaged pulse')
        plt.plot(pulse_list[diff_list.argmax()],label='Worst outlayer in current process')
        plt.title("Pulse averaging process")
        plt.legend(loc="best")
        plt.show(block=False)
        answer = input("Suppress 10% worst [y/n]?")
        if answer=='y':
            pulse_list = pulse_list[(diff_list<np.percentile(diff_list,90))]
        else:
            satisfied=True
            
    return mean_pulse

# ############################################################
# Function to estimate a noise spectrum average
# ############################################################
def accumulate_noise_spectra(noise_records,pulse_length,abs_mean=False,rebin=1,normalize=False,dt=6.4e-6):
    '''Accumulates noise spectra from pulse free data streams
    
    Arguments:
        - noise_records: input pulse-free records
        - pulse_length: length of the records, fft to perform
        - abs_mean: option to average the abs values instead of the squares
        - rebin: rebinning factor to apply on the final spectrum
        - normalize: option to normalize the spectrum in proper A/rHz units
        - dt: sampling rate of the data (only relevant in case of normalize option)
        
    Returns: average noise spectra in a np vector of length pulse_length 
    '''
    noise_spectrum_tot = np.zeros(int(pulse_length/rebin))
    nb_noise_estimates=0
    for stream_data in noise_records:
        for i in range(int(len(stream_data)/pulse_length)):
            if abs_mean:
                noise_spectrum_tot+=(abs(np.fft.fft(stream_data[i][:pulse_length]))).reshape(-1,rebin).mean(1)
            else:
                noise_spectrum_tot+=(abs(np.fft.fft(stream_data[i][:pulse_length]))**2).reshape(-1,rebin).mean(1)
            nb_noise_estimates+=1
    print("Number of records used in noise spectrum calibration:",nb_noise_estimates)
    noise_spectrum_tot/=nb_noise_estimates
    if not abs_mean:
        noise_spectrum_tot = np.sqrt(noise_spectrum_tot)
    if normalize:
        noise_spectrum_tot*=np.sqrt(2*dt/pulse_length)
    return noise_spectrum_tot

# ############################################################
# Function to estimate a noise spectrum average
# ############################################################
def accumulate_noise_spectra_ini(record_data,pulse_length,abs_mean=False,rebin=1,normalize=False,dt=6.4e-6):
    '''Accumulates noise spectra from pulse free data streams
    
    Arguments:
        - record_data: input pulse-free data stream
        - pulse_length: length of the fft to perform
        - abs_mean: option to average the abs values instead of the squares
        - rebin: rebinning factor to apply on the final spectrum
        - normalize: option to normalize the spectrum in proper A/rHz units
        - dt: sampling rate of the data (only relevant in case of normalize option)
        
    Returns: average noise spectra in a np vector of length pulse_length 
    '''
    noise_spectrum_tot = np.zeros(int(pulse_length/rebin))
    nb_noise_estimates=0
    for stream_data in record_data:
        for i in range(int(len(stream_data)/pulse_length)):
            if abs_mean:
                noise_spectrum_tot+=(abs(np.fft.fft(stream_data[i*pulse_length:(i+1)*pulse_length]))).reshape(-1,rebin).mean(1)
            else:
                noise_spectrum_tot+=(abs(np.fft.fft(stream_data[i*pulse_length:(i+1)*pulse_length]))**2).reshape(-1,rebin).mean(1)
            nb_noise_estimates+=1
    print("Number of segments used in noise spectrum calibration:",nb_noise_estimates)
    noise_spectrum_tot/=nb_noise_estimates
    if not abs_mean:
        noise_spectrum_tot = np.sqrt(noise_spectrum_tot)
    if normalize:
        noise_spectrum_tot*=np.sqrt(2*dt/pulse_length)
    return noise_spectrum_tot

# ############################################################
# Function to compute an optimal filter from a pulse template and a noise spectrum
# ############################################################
def compute_optimal_filter(pulse_template,noise_spectrum,energy):
    """Function to compute an optimal filter from a pulse template and a noise spectrum.
    Optimal filter is normalized to return energy when dot producted with the template
    
    Arguments:
        - pulse_template: pulse template to use in optimal filter generation
        - noise_spectrum: noise spectrum to use in optimal filter generation
        - energy: energy of the pulse template 
    """
    pulse_spectrum = np.fft.fft(pulse_template)
    time_optimal_filter = pulse_spectrum/noise_spectrum**2
    time_optimal_filter[0]=0
    cutted_time_filter=np.real(np.fft.ifft(time_optimal_filter))
    cutted_time_filter/=np.dot(pulse_template,cutted_time_filter)/energy
    
    return cutted_time_filter


# ############################################################
# Function to perform the "jitter parabola fit"
# ############################################################
def do_pulse_jitter(opt_filter,stream_data,pulse_time):
    '''Performs the +-1 pulse jitter parabola technique and returns both the maximum and its phase and
    the potentially corrected pulse time
    
    Arguments:
        - opt_filter: optimal filter
        - stream_data: data_stream to analyze
        - pulse_time: guess arrival time (if at some )
    '''
    phase_offset = 0 
    pulse_length = len(opt_filter)
    while True:
        energy1 = np.dot(opt_filter,stream_data[pulse_time-1:pulse_time-1+pulse_length])
        energy2 = np.dot(opt_filter,stream_data[pulse_time:pulse_time+pulse_length])
        energy3 = np.dot(opt_filter,stream_data[pulse_time+1:pulse_time+1+pulse_length])
        A = 0.5*(energy1-2*energy2+energy3)
        B = 0.5*(energy3-energy1)
        C = energy2
        energy = C-.25*B**2/A
        if energy1>energy2:
            pulse_time-=1
            phase_offset-=1
        elif energy3>energy2:
            pulse_time+=1
            phase_offset+=1
        else:
            pulse_phase = -.5*B/A
            break
    return energy,pulse_phase+phase_offset,pulse_time


# ############################################################
# Gaussian function
# ############################################################
def gauss(x,*p):
    '''Gaussian function for fits
    '''
    A,mu,sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


# ############################################################
# Generic function to fit a histogram
# ############################################################
def hist_and_fit(array_to_fit,bins,show=False,log=False,xlabel="",fit_function=gauss,figfile=None,c=(71./255,144./255,195./255)):
    '''Plots an histogram and fit it with a gaussian
    
    Arguments:
        - array_to_fit: array containing the data to study
        - bins: bins parameter of np.histogram
        - log: option to plot in log scale
        - xlabel: label of x axis for plot
        - fit_function: callable for the fit function
        - figfile: if given, path to file to save the fit plot
        - c: RGB face color of the histogram
    '''
    hist,bins = np.histogram(array_to_fit,bins=bins)
    bin_centers = (bins[:-1] + bins[1:])/2
    p0=[len(array_to_fit)/(np.sqrt(2*np.pi)*array_to_fit.std()),array_to_fit.mean(),array_to_fit.std()]
    coeff = curve_fit(fit_function, bin_centers, hist, p0=p0,maxfev=10000000)[0]
    if show or figfile is not None:
        axe_fit = np.arange(bins[0],bins[-1],0.01*(bins[-1]-bins[0]))
        hist_fit = fit_function(axe_fit, *coeff)
        plt.hist(array_to_fit,bins=bins,histtype='stepfilled',facecolor=c)
        plt.plot(axe_fit,hist_fit,c='r',linewidth=3)
        if log:
            plt.yscale('log')
        plt.xlabel(xlabel)
        plt.ylabel("Occurences")
        print("Fitted function with parameters:",coeff)
        if fit_function==gauss:
            print("  - mu: {0:5.3f}".format(coeff[1]))
            print("  - sigma: {0:5.3f} (FWHM = {1:5.3f})".format(coeff[2],2.355*coeff[2]))
            print("  - A: {0:5.3f}".format(coeff[0]))
        if figfile is not None:
            plt.savefig(figfile,bbox_inches='tight')
            print("Figure saved to"+figfile)
        if show:
            plt.show()
        plt.clf()
        
    return coeff

# ############################################################
# Function for phase correction
# ############################################################
def apply_phase_correction(energies,phase,phase_correction,verbose=True):
    '''Applies a phase correction (should be merged with baseline correction)
    
    Arguments:
        - energies: input energies
        - phase: phase data
        - phase_correction: order of the phase correction polynom
        - verbose: option to show the corrected values
    '''
    phase_correc_poly_coeff = np.polyfit(phase,energies,phase_correction)
    phase_correc_poly = np.poly1d(phase_correc_poly_coeff)
    plt.scatter(phase,energies)
    plt.plot(np.sort(phase),phase_correc_poly(np.sort(phase)),c='r',linewidth=2)
    plt.xlabel('Phase')
    plt.ylabel('Energy [keV]')
    plt.show()
    
    #Correct for correlation
    corrected_energies=energies-phase_correc_poly(phase)+phase_correc_poly(phase.mean())
    if verbose:
        plt.scatter(phase,corrected_energies)
        plt.xlabel('Phase')
        plt.ylabel('Energy [keV]')
        plt.show()
    
    return corrected_energies

# ############################################################
# Function for baseline correction
# ############################################################
def apply_baseline_correction(energies,baseline,baseline_correction,verbose=True):
    '''Applies a phase correction
    
    Arguments:
        - energies: input energies
        - baseline: baseline data
        - baseline_correction: order of the baseline correction polynom
        - verbose: option to show the corrected values
    '''
    #Fit baseline correction
    baseline_correc_poly_coeff = np.polyfit(baseline,energies,baseline_correction)
    baseline_correc_poly = np.poly1d(baseline_correc_poly_coeff)
    plt.scatter(baseline,energies)
    plt.plot(np.sort(baseline),baseline_correc_poly(np.sort(baseline)),c='r',linewidth=2)
    plt.xlim((min(baseline)*0.99,1.01*max(baseline)))
    plt.xlabel('Baseline')
    plt.ylabel('Energy [keV]')
    plt.show()
    
    #Correct for correlation
    energies-=baseline_correc_poly(baseline)-baseline_correc_poly(baseline.mean())
    if verbose:
        plt.scatter(baseline,energies)
        plt.xlim((min(baseline)*0.99,1.01*max(baseline)))
        plt.xlabel('Baseline')
        plt.ylabel('Energy [keV]')
        plt.show()
    
    return energies


# ############################################################
# Function to Estimate threshold level from initial noise fluctuations
# ############################################################
def get_threshold(noise_record,nsigmas=7,verbose=False,forced_threshold=None):
    """Define the threshold level from noise data
    
    Arguments:
        - noise_record: array of stream values to be triggered
        - nsigmas: significance of the pulse detection with respect to the initial stream fluctuations
        - verbose: option to print and plot additional information, not relevant for nominal script use
        - forced_threshold: if given, use this value as threshold instead of computing it from the start of the data
    """
    if forced_threshold is not None:
        threshold = forced_threshold
    else:
        noise = noise_record
        derivated_noise = np.convolve(noise,[1,-1],mode="valid")
        threshold = nsigmas*np.std(derivated_noise)
        if verbose:
            print("Chosen threshold:",threshold)
        if verbose:
            plt.plot(derivated_noise)
            plt.title("Data used for threshold estimation")
            plt.show()
    return threshold
        

# ############################################################
# Function to trigger pulses in a stream
# ############################################################
def trigger_pulses(event_records,threshold,pulse_length=2048,prebuffer=PREBUFF,
                   verbose=False,jitter_margin=0,silent=False):
    """Detects individual pulses in a stream and returns the list of pulses
    
    Arguments:
        - event_records: array of event records to be triggered
        - threshold: threshold level for pulse detection
        - pulse_length: length of the pulses to detect
        - prebuffer: length of data to leave before each pulse
        - verbose: option to print and plot additional information, not relevant for nominal script use
        - jitter_margin: number of points to add on each side of triggered pulses to allow jitter correction
    """
    if not silent:
        print("WARNING: this function was written without any care for pileup detection")

    pulses=np.array((len(event_records, pulse_length)))
    pulse_time=np.array(len(event_records))

    for i in range(len(event_records)):
        record=event_records[i]
        # Derivate data
        derivated_record = np.convolve(record,[1,-1],mode="valid")
        if verbose:
            plt.plot(derivated_record)
            plt.plot(threshold*np.ones_like(derivated_record))
            plt.title("Comparison of pulses and threshold")
            plt.show()

        i_threshold=np.where(derivated_record>threshold)[0]
        pulse_time[i] = i_threshold-prebuffer-1
        if verbose:
            print('Pulse detected at time',pulse_time[i])
        pulses[i]=record[pulse_time-jitter_margin:pulse_time+pulse_length+jitter_margin]
            
    # Overplot triggered pulses if requested 
    if verbose:
        plt.plot(stream_data)
        for pulse,pulse_time in zip(pulses,pulse_times):
            plt.plot(np.arange(pulse_time-jitter_margin,pulse_time+jitter_margin+pulse_length,1),pulse,c='green')
        plt.show()
        
        for pulse in pulses:
            plt.plot(pulse)
        plt.show()
    
    # Print trigger statistics and return
    if not silent:
        print("Triggered {0:6d} pulses".format(len(pulse_times)))
    return pulses,pulse_times,threshold


# ############################################################
# Function to perform energy reconstruction
# ############################################################
def energy_reconstruction(pulse_list,optimal_filter,pulse_time,prebuffer=PREBUFF,prebuff_exclusion=10,
                          verbose=False):
    """Perform energy reconstruction on list of pulses, including jitter correction.
    Also compute baseline value before each pulse
    
    Arguments:
        - pulse_list: list of pulses to reconstruct
        - optimal_filter: filter to use for reconstruction
        - prebuffer: length of prebuffer data available before each pulse
        - prebuff_exclusion: number of points to remove from the buffer in baseline estimation 
        
    Returns: pulse_template 
        - pulse_template: template created from the average of detected pulses
    """
    # Iterate over available pulses
    energies = []
    phases = []
    times = []
    baselines = []
    for pulse in pulse_list:
        # Perform reconstruction with jitter
        e,p,t = do_pulse_jitter(optimal_filter, pulse, pulse_time)
        energies.append(e)
        phases.append(p)
        times.append(t)
        
        # Estimate baseline
        baselines.append(pulse[:prebuffer-prebuff_exclusion].mean())
        
        
    energies = np.array(energies)
    phases = np.array(phases)
    times = np.array(times)
    baselines = np.array(baselines)
    
    # Return energies, phases and times
    return energies, phases, times, baselines
            

if __name__ == '__main__':
    # ############################################################
    # Parse arguments
    # ############################################################
    
    parser = argparse.ArgumentParser(description='Script to analyze data generated by the DRE DEMUX DM')
    parser.add_argument("--threshold", type=float, help="Level of threshold for pulse detection (in sigmas)", default=7)
    parser.add_argument("--pulse_length", "-p", type=float, help="Pulse length for reconstruction", default=2048)
    parser.add_argument("--ignore_outlayers", "-i", help="Option to ignore outlayers in pulse template generation", action="store_true")
    parser.add_argument("--tes_noise", "-t", help="Option to add the TES noise to the DAC noise spectrum", action="store_true")
    parser.add_argument("--bcorr", type=float, help="Order of polynomial baseline correction", default=3)
    parser.add_argument("--pcorr", type=float, help="Order of polynomial arrival phase correction", default=8)
    parser.add_argument("--show", "-s", help="Option to show diagnostic plots", action="store_true")
    parser.add_argument("--verbose", help="Option to print additional information", action="store_true")
    parser.add_argument("--directory", "-d", type=str, help="PATH to directory with data", default=DRE_DATA_ROOT)
    parser.add_argument("--outdir", "-o", type=str, help="PATH to directory where plots should be saved", default=None)
    parser.add_argument("--ep_params", type=str, help="PATH to directory where DRE EP parameters were saved", default=None)
    args = parser.parse_args()
    
    do_plots = args.show or (args.outdir is not None)
    

    print("\nUsing data in directory",args.directory)
    # ############################################################
    # Pulse template calibration
    # ############################################################
    
    print("\nPerforming pulse template calibration...")
    
    # Load pulse data to be used for template calibration
    pulse_data=-get_records(args.directory+"event_records.dat", verbose=args.verbose)
    #JUST FOR DEBUG
    #pulse_data=pulse_data[:int(len(pulse_data)*0.7)]

    print("===>> ", len(pulse_data))
    print("===>> ", len(pulse_data[10]))
    plt.plot(pulse_data[23])
    plt.show()

    # Define threshold level
    threshold=get_threshold(noise_record,nsigmas=7,verbose=args.verbose,forced_threshold=None)

    # Trigger individual pulses
    pulse_list,pulse_times=trigger_pulses(event_records,threshold,pulse_length=2048,prebuffer=PREBUFF,
                   verbose=False,jitter_margin=0,silent=False)

    #pulse_list,pulse_times,threshold = trigger_pulses(pulse_data,args.threshold,threshold_estimate=THRESHOLD_ESTIMATE,verbose=args.verbose,
    #                                                  pulse_length=args.pulse_length,prebuffer=PREBUFF)
    
    # Generate pulse template as average of detected pulses
    pulse_template = pulse_average(pulse_list,ignore_outlayers=args.ignore_outlayers)
    if args.verbose:
        print("Pulse template:",pulse_template)
    if do_plots:
        plt.plot((np.arange(args.pulse_length)-200)*6.4e-3,pulse_template)
        plt.title('Pulse template')
        plt.xlabel('Time [ms]')
        plt.ylabel("-ADU")
        if args.outdir is not None:
            plt.savefig(args.outdir+"/pulse_template.pdf",bbox_inches='tight')
        if args.show:
            plt.show()
        plt.clf()
    
    # ############################################################
    # Comparison with xifusim data
    # ############################################################
    
    print("\nComparing xifusim and dre data...")
    
    # Load xifusim template
    xifusim_time,xifusim_template = np.load(XIFUSIM_TEMPLATE)
    
    # Determine scaling factor between xifusim and DRE units with the baseline value
    dre_template = -pulse_template
    baseline_scaling = dre_template[0]/xifusim_template[0]
    PH_scaling = (dre_template[0]-dre_template.min())/(xifusim_template[0]-xifusim_template.min())
    print("Difference between average baseline and PH scaling: {0:5.3f}%".format((PH_scaling/baseline_scaling-1)*100))
    
    if do_plots:
        # Plot xifusim pulse vs. variability of individual pulses
        for pulse in pulse_list[:20]:
            plt.plot((np.arange(args.pulse_length)-201.5)*6.4e-3,-pulse)
        plt.plot(xifusim_time*1e3,xifusim_template*baseline_scaling,label="xifusim",c="red",linewidth=3)
        plt.title('DRE pulses vs. xifusim pulse template')
        plt.xlabel('Time [ms]')
        plt.ylabel("ADU")
        plt.legend(loc="best")
        if args.outdir is not None:
            plt.savefig(args.outdir+"/pulse_template_comp_xifusim.pdf",bbox_inches='tight')
        if args.show:
            plt.show()
        plt.clf()
    
        # Show difference with respect to overall pulses
        xifusim_template_decimated = np.convolve(xifusim_template,1/128.*np.ones(128),mode="valid")[79::128][:1000] # Fitted best phase by hand
        plt.plot(xifusim_template_decimated*baseline_scaling,label='xifusim')
        plt.plot(dre_template[195:1195],label='dre')
        plt.plot(xifusim_template_decimated*baseline_scaling - dre_template[195:1000+195],label='difference')
        plt.title('Direct difference between both templates')
        plt.legend(loc="best")
        if args.outdir is not None:
            plt.savefig(args.outdir+"/pulse_template_abs_diff_xifusim.pdf",bbox_inches='tight')
        if args.show:
            plt.show()
        plt.clf()
        
        # Relative difference
        plt.plot((xifusim_template_decimated*baseline_scaling - dre_template[195:1000+195])/dre_template[195:1000+195]*100,label='difference')
        plt.ylabel('Relative difference [%]')
        plt.title('Relative difference between both templates')
        if args.outdir is not None:
            plt.savefig(args.outdir+"/pulse_template_relative_diff_xifusim.pdf",bbox_inches='tight')
        if args.show:
            plt.show()
        plt.clf()
        
        # Compare power spectra
        xifusim_PS = abs(np.fft.rfft(xifusim_template_decimated*baseline_scaling))**2
        DRE_PS = abs(np.fft.rfft(dre_template[195:1195]))**2
        PS_freq = np.fft.fftfreq(1000,6.4e-6)[:501]
        plt.loglog(PS_freq,xifusim_PS,label="xifusim")
        plt.loglog(PS_freq,DRE_PS,label="dre")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("PSD [AU]")
        plt.title('Comparison of power spectra')
        plt.legend(loc="best")
        if args.outdir is not None:
            plt.savefig(args.outdir+"/pulse_template_spectrum_diff_xifusim.pdf",bbox_inches='tight')
        if args.show:
            plt.show()
        plt.clf()
        
        
    # ############################################################
    # Noise spectrum calibration
    # ############################################################
    
    print("\nPerforming noise spectrum calibration...")
    
    # Load pulse-free data to generate noise spectrum
    noise_data = np.load(args.directory+"noise.npy")[1]
    if args.verbose:
        plt.plot(noise_data)
        plt.xlabel("Sample")
        plt.ylabel("ADU")
        plt.show()
        
    # Compute average noise spectrum
    noise_spectrum = accumulate_noise_spectra([noise_data[:args.pulse_length*75]], args.pulse_length, normalize=True)
    frequencies = np.fft.fftfreq(args.pulse_length,6.4e-6)[1:int(args.pulse_length/2)]
    if args.verbose:
        print("Noise spectrum:",noise_spectrum)
    if args.show:
        plt.loglog(frequencies,noise_spectrum[1:int(args.pulse_length/2)])
        plt.title('Average noise spectrum')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel("ADU/rHz")
        plt.show()
            
            
    # ############################################################
    # Comparison with TES noise
    # ############################################################
    
    print("\nComparing TES noise and DRE noise...")
    
    # Load TES noise spectrum
    tes_noise = fits.getdata(TES_NOISE_SPECTRUM,"SPEC{0:04d}".format(args.pulse_length))["CSD"]*baseline_scaling
    total_noise = np.sqrt(tes_noise**2+noise_spectrum**2)
    
    # Compare DAC noise with TES noise
    if do_plots:
        plt.loglog(frequencies,tes_noise[1:int(args.pulse_length/2)],label="TES noise")
        plt.loglog(frequencies,noise_spectrum[1:int(args.pulse_length/2)],label="DAC noise")
        plt.loglog(frequencies,total_noise[1:int(args.pulse_length/2)],label="Total noise")
        plt.title('Average noise spectrum')
        plt.title("DAC noise vs TES noise")
        plt.legend(loc="best")
        if args.outdir is not None:
            plt.savefig(args.outdir+"/noise_spectrum_vs_TES.pdf",bbox_inches='tight')
        if args.show:
            plt.show()
        plt.clf()
    
    # ############################################################
    # Compute optimal filter
    # ############################################################
    
    print('\nComputing optimal filter...')
    optimal_filter = compute_optimal_filter(pulse_template,noise_spectrum,energy=7.)
    print("Check optimal filter normalization: {0:5.3f}keV".format(np.dot(optimal_filter,pulse_template)))
    
    
    # ############################################################
    # Pulse reconstruction
    # ############################################################
    
    print('\nReconstructing pulses...')
    NONLINEAR_FACTOR=get_nonlinear_factor(PIXEL_DATA_DIR, verbose=args.verbose)

    # Load pulse data containing the events to reconstruct
    pulse_data = -np.load(args.directory+"measure.npy")[1]
    
    # Trigger individual pulses
    pulse_list,pulse_times,_ = trigger_pulses(pulse_data,args.threshold,threshold_estimate=THRESHOLD_ESTIMATE,prebuffer=PREBUFF,verbose=args.verbose,
                                            pulse_length=args.pulse_length,forced_threshold=threshold,jitter_margin=JITTER_MARGIN)
    
    # Compute first energy estimates
    energies,phases,_,baselines = energy_reconstruction(pulse_list, optimal_filter, JITTER_MARGIN)
    coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, xlabel="Error before any correction [eV]")
    print("Resolution without TES noise (prior correction): {0:5.3f}eV".format(coeffs[2]*2.355*NONLINEAR_FACTOR))
        
    # Apply baseline correction
    energies = apply_baseline_correction(energies, baselines, args.bcorr, verbose=args.verbose)
    coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, xlabel="Error after baseline correction [eV]")
    print("Resolution without TES noise (after baseline correction): {0:5.3f}eV".format(coeffs[2]*2.355*NONLINEAR_FACTOR))
    
    # Apply phase correction
    energies = apply_phase_correction(energies, phases, args.pcorr, verbose=args.verbose)
    if args.outdir is not None:
        coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, figfile=args.outdir+"/e_error_noTESnoise.pdf", xlabel="Final reconstruction error [eV]")
    else:
        coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, xlabel="Final reconstruction error [eV]")
    eres_notesnoise = coeffs[2]*2.355*NONLINEAR_FACTOR
    print("Final resolution without TES noise: {0:5.3f}+-{1:5.3f}eV".format(eres_notesnoise,eres_notesnoise/(np.sqrt(2.*len(energies)))))
    
    # Compare error from Paul and mine
    if (args.show or args.outdir is not None) and os.path.exists(args.directory+"events.npy"):
        events = np.load(args.directory+"events.npy")[-1]
        if len(events)!=len(energies):
            print("\nWARNING! Did not detect the same number of pulses as Paul's code -> skip comparison")
        else:
            plt.scatter(events-7000,energies*1000-7000)
            plt.xlabel("Energy reconstruction error by IRAP [eV]")
            plt.ylabel("Energy reconstruction error by Philippe [eV]")
            if args.outdir is not None:
                plt.savefig(args.outdir+"/pprec_vs_IRAP.pdf",bbox_inches='tight')
            if args.show:
                plt.show()
            plt.clf()

    # ############################################################
    # Comparison using DRE template
    # ############################################################
    if args.ep_params is not None:
        if not os.path.exists(args.directory+"events.npy"):
            raise RuntimeError("Tried to apply IRAP optimal filter with no corresponding data -> What's the point?")
        opt_filt_irap = np.genfromtxt(args.ep_params+"/Pattern.txt")
        b_irap,events = np.load(args.directory+"events.npy")[1:]
        e_irap,p_irap,_,_ = energy_reconstruction(pulse_list, opt_filt_irap, JITTER_MARGIN)
        e_irap*=7000/e_irap.mean()
        e_irap = apply_baseline_correction(e_irap, b_irap, args.bcorr, verbose=args.verbose)
        e_irap = apply_phase_correction(e_irap, p_irap, args.bcorr, verbose=args.verbose)
        if args.show or args.outdir is not None:
            plt.scatter(events-7000,e_irap-7000)
            plt.xlabel("Energy reconstruction error by IRAP [eV]")
            plt.ylabel("Energy reconstruction error by Philippe with IRAP template [eV]")
            if args.outdir is not None:
                plt.savefig(args.outdir+"/pprec_vs_IRAP_IRAPtemplate.pdf",bbox_inches='tight')
            if args.show:
                plt.show()
            plt.clf()
        
        if args.outdir is not None:
            coeffs = hist_and_fit(e_irap-7000, 100, show=args.show, figfile=args.outdir+"/e_error_noTESnoise_IRAPtemplate.pdf", xlabel="Final reconstruction error with IRAP template [eV]")
        else:
            coeffs = hist_and_fit(e_irap-7000, 100, show=args.show, xlabel="Final reconstruction error with IRAP template [eV]")
        eres_irap = coeffs[2]*2.355*NONLINEAR_FACTOR
        print("Final resolution with IRAP template: {0:5.3f}+-{1:5.3f}eV".format(eres_irap,eres_irap/(np.sqrt(2.*len(energies)))))

    
    # ############################################################
    # Pulse reconstruction with OF also containing TES noise
    # ############################################################
    
    print("\nReconstruction with OF including TES noise...")
    
    # Compute optimal filter, now including TES noise
    optimal_filter_tot = compute_optimal_filter(pulse_template,total_noise,7.)
    print("Check optimal filter normalization: {0:5.3f}".format(np.dot(optimal_filter_tot,pulse_template)))
    
    # Compute first energy estimates
    energies,phases,_,baselines = energy_reconstruction(pulse_list, optimal_filter_tot, JITTER_MARGIN)
    coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, xlabel="Error before any correction [eV]")
    print("Resolution with TES noise (prior correction): {0:5.3f}ev".format(coeffs[2]*2.355*NONLINEAR_FACTOR))
        
    # Apply baseline correction
    energies = apply_baseline_correction(energies, baselines, args.bcorr, verbose=args.verbose)
    coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, xlabel="Error after baseline correction [eV]")
    print("Resolution with TES noise (after baseline correction): {0:5.3f}eV".format(coeffs[2]*2.355*NONLINEAR_FACTOR))
    
    # Apply phase correction
    energies = apply_phase_correction(energies, phases, args.pcorr, verbose=args.verbose)
    if args.outdir is not None:
        coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, figfile=args.outdir+"/e_error_TESnoise.pdf", xlabel="Final reconstruction error [eV]")
    else:
        coeffs = hist_and_fit((energies-7)*1000, 100, show=args.show, xlabel="Final reconstruction error [eV]")
    eres_tesnoise = coeffs[2]*2.355*NONLINEAR_FACTOR
    print("Final resolution with TES noise: {0:5.3f}+-{1:5.3f}eV".format(eres_tesnoise,eres_tesnoise/(np.sqrt(2.*len(energies)))))
        
    # Compare error from Paul and mine
    if False and args.show and os.path.exists(args.directory+"events.npy"):
        events = np.load(args.directory+"events.npy")[-1]
        plt.scatter(events-7000,energies*1000-7000)
        plt.xlabel("Energy reconstruction error by Paul [eV]")
        plt.ylabel("Energy reconstruction error by Philippe [eV]")
        plt.show()
        
        
    
