'''
Created on 26 juil. 2019
@author: philippepeille

Updated on 3 Sept. 2019 by Laurent Ravera to match the DRE processing environnement
'''
import os
import numpy as np
import argparse
import general_tools
import matplotlib.pyplot as plt
import general_tools
from astropy.io import fits
from scipy.optimize.minpack import curve_fit

PREBUFF=180
JITTER_MARGIN=10


# ############################################################
# Function to read data records from fits files
# ############################################################
def get_records_from_fits(filename, verbose=False):
    '''Reads sample records from a fits file

    Arguments:
        - filename: name of the file containing the sample records
        - verbose: option to print the nonlinear factor

    Returns: an array containing the records
    '''

    hdul = fits.open(filename)
    data=hdul[1].data
    t=data['Timestamp']
    chid=data['channelNum']
    pixid=data['pixelNum']
    i=data['i']
    q=data['q']
    module=np.sqrt(i.astype("float")**2+q.astype("float")**2)
    if verbose:
        print("  Informations of FITS file:")
        print("    Date:    ", hdul[1].header['DATE'])
        print("    Origin:  ", hdul[1].header['ORIGIN'])
        print("    Project: ", hdul[1].header['INSTRUME'])
        print("    Number of records: ", len(i))
        print("    Length of records: ", int(len(i[0])))
    
    return(chid, pixid, module, t)
    
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
    return NLF

# ############################################################
# Function to create pulse average
# ############################################################
def pulse_average(pulse_list,ignore_outlayers=False):
    """Creates pulse template from a set of pulses. Outlayers can be manually rejected.
    
    Arguments:
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
        plt.legend(loc="best", prop=dict(size=7))
        plt.show(block=False)
        answer = input("  Suppress 10% worst [y/n]?")
        if answer=='y':
            pulse_list = pulse_list[(diff_list<np.percentile(diff_list,90))]
        else:
            satisfied=True
            
    return mean_pulse


# ############################################################
# Function to estimate a noise spectrum average
# ############################################################
def accumulate_noise_spectra(noise_records,abs_mean=False,rebin=1,normalize=False,dt=6.4e-6):
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
    nb_records=len(noise_records)
    pulse_length=len(noise_records[0])
    noise_spectrum_tot = np.zeros(int(pulse_length/rebin))
    nb_noise_estimates=0
    for i_record in range(nb_records):
        record=noise_records[i_record]
        if abs_mean:
            noise_spectrum_tot+=(abs(np.fft.fft(record))).reshape(-1,rebin).mean(1)
        else:
            noise_spectrum_tot+=(abs(np.fft.fft(record))**2).reshape(-1,rebin).mean(1)
        nb_noise_estimates+=1
    print("  Number of records used in noise spectrum calibration:",nb_noise_estimates)
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
def do_pulse_jitter(opt_filter,pulse_record):
    '''Performs the +-1 pulse jitter parabola technique and returns both the maximum and its phase and
    the potentially corrected pulse time
    
    Arguments:
        - opt_filter: optimal filter
        - pulse_record: data_stream to analyze
    '''
    phase_offset = 5 
    pulse_length = len(opt_filter)
    while True:
        if phase_offset<0 or phase_offset+pulse_length>len(pulse_record):
            print('Problem to find the phase of the pulse!')
            energy=-1
            break
        energy1 = np.dot(opt_filter,pulse_record[phase_offset:phase_offset+pulse_length])
        energy2 = np.dot(opt_filter,pulse_record[phase_offset+1:phase_offset+pulse_length+1])
        energy3 = np.dot(opt_filter,pulse_record[phase_offset+2:phase_offset+pulse_length+2])
        A = 0.5*(energy1-2*energy2+energy3)
        B = 0.5*(energy3-energy1)
        C = energy2
        energy = C-.25*B**2/A
        pulse_phase = -.5*B/A

        if energy1>energy2:
            #pulse_time-=1
            phase_offset-=1
        elif energy3>energy2:
            #pulse_time+=1
            phase_offset+=1
        else:
            pulse_phase = -.5*B/A
            break

    return energy,pulse_phase+phase_offset


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
def hist_and_fit(array_to_fit,bins,fit_function=gauss,c=(71./255,144./255,195./255)):
    '''Plots an histogram and fit it with a gaussian
    
    Arguments:
        - array_to_fit: array containing the data to study
        - bins: bins parameter of np.histogram
        - xlabel: label of x axis for plot
        - fit_function: callable for the fit function
        - c: RGB face color of the histogram
    '''
    hist,bins = np.histogram(array_to_fit,bins=bins)
    bin_centers = (bins[:-1] + bins[1:])/2
    p0=[len(array_to_fit)/(np.sqrt(2*np.pi)*array_to_fit.std()),array_to_fit.mean(),array_to_fit.std()]
    coeff = curve_fit(fit_function, bin_centers, hist, p0=p0,maxfev=10000000)[0]

    axe_fit = np.arange(bins[0],bins[-1],0.01*(bins[-1]-bins[0]))
    hist_fit = fit_function(axe_fit, *coeff)

    print("Fitted function with parameters:",coeff)
    if fit_function==gauss:
        print("  - mu: {0:5.3f}".format(coeff[1]))
        print("  - sigma: {0:5.3f} (FWHM = {1:5.3f})".format(coeff[2],2.355*coeff[2]))
        print("  - A: {0:5.3f}".format(coeff[0]))
        
    return coeff, array_to_fit, bins, axe_fit, hist_fit


# ############################################################
# Function for phase correction
# ############################################################
def apply_phase_correction(energies,phase,phase_correction):
    '''Applies a phase correction (should be merged with baseline correction)
    
    Arguments:
        - energies: input energies
        - phase: phase data
        - phase_correction: order of the phase correction polynom
        - fig: figure id to include the phase correction plots
        - i_subplot: index of the sub_plot in the figure
    '''
    phase_correc_poly_coeff = np.polyfit(phase,energies,phase_correction)
    phase_correc_poly = np.poly1d(phase_correc_poly_coeff)
    corrected_energies=energies.copy()-phase_correc_poly(phase)+phase_correc_poly(phase.mean())
   
    return corrected_energies, phase_correc_poly

# ############################################################
# Function for baseline correction
# ############################################################
def apply_baseline_correction(energies,baseline,baseline_correction):
    '''Applies a phase correction
    
    Arguments:
        - energies: input energies
        - baseline: baseline data
        - baseline_correction: order of the baseline correction polynom
        - fig: figure id to include the phase correction plots
        - i_subplot: index of the sub_plot in the figure
    '''
    #Fit baseline correction
    baseline_correc_poly_coeff = np.polyfit(baseline,energies,baseline_correction)
    baseline_correc_poly = np.poly1d(baseline_correc_poly_coeff)
    energies_corrected=energies.copy()-baseline_correc_poly(baseline)+baseline_correc_poly(baseline.mean())
    
    return energies_corrected, baseline_correc_poly


# ############################################################
# Function to perform energy reconstruction
# ############################################################
def energy_reconstruction(pulse_list,optimal_filter,prebuffer=PREBUFF,prebuff_exclusion=10,
                          verbose=False):
    """Perform energy reconstruction on list of pulses, including jitter correction.
    Also compute baseline value before each pulse
    
    Arguments:
        - pulse_list: list of pulses to reconstruct
        - optimal_filter: filter to use for reconstruction
        - prebuffer: length of prebuffer data available before each pulse
        - prebuff_exclusion: number of points to remove from the buffer in baseline estimation 
        
    Returns: pulse_template 
        - 
    """
    # Iterate over available pulses
    energies = []
    phases = []
    times = []
    baselines = []
    for i_pulse in range (len(pulse_list)):
        # Perform reconstruction with jitter
        e,p = do_pulse_jitter(optimal_filter, pulse_list[i_pulse])
        if e!=-1:  # an optimal phase has been found
            energies.append(e)
            phases.append(p)
            # Estimate baseline
            baselines.append(pulse_list[i_pulse][:prebuffer-prebuff_exclusion].mean())
             
    energies = np.array(energies)
    phases = np.array(phases)
    baselines = np.array(baselines)
    
    # Return energies, phases and times
    return energies, phases, baselines


# ############################################################
# Function to perform all the operations needed to compute the optimal filters
# ############################################################
def do_EP_filter(file_noise, file_pulses, file_xifusim_template, file_xifusim_tes_noise, plotdirname, verbose=False, do_plots=True):
    """Perform the operations to compute the optimal filters (with and without tes noise)
    and compares the results with xifusim.
    
    Arguments:
        - file_noise: fits file containing DRE noise data
        - file_pulses: fits file containing DRE pulses data
        - file_xifusim_template: file containing xifusim template
        - file_xifusim_tes_noise: file containing tes noise compute by xifusim
        - plotdirname: location of plotfiles
        - verbose: if True some informations are printed (Default=False)
        - do_plots: if True a plot is done with all the intermediate results (Default=True)
        
    Returns: optimal_filter (no TES noise), optimal_filter_tot (with TES noise)       
    """
 
    # ############################################################
    # Noise spectrum calibration
    # ############################################################
    print("\nPerforming noise spectrum calibration...")

    # Load pulse-free data to generate noise spectrum
    print("  Loading noise data from file ", file_noise)
    _, pix, module_n, _ = get_records_from_fits(file_noise, verbose=verbose)
    i_pixtest = np.where(pix==40)[0] # keeping only records from test pixel
    noise_data = module_n[i_pixtest] 
    record_length=len(noise_data[0])
    print("  Record length = {0:4d}".format(record_length))
        
    # Compute average noise spectrum
    noise_spectrum = accumulate_noise_spectra(noise_data, normalize=True)
    np.save('noise_spectrum.npy', noise_spectrum)
    frequencies = np.fft.fftfreq(record_length,6.4e-6)[1:int(record_length/2)]
    if verbose:
        print("  Noise spectrum:",noise_spectrum)

    # ############################################################
    # Pulse template calibration
    # ############################################################
        print("\nPerforming pulse template calibration...")
    
    # Load pulse data to be used for template calibration
    print("  Loading calibration pulse data from file ", file_pulses)
    _, _, pulse_list, pulse_times=get_records_from_fits(file_pulses, verbose=verbose)
    print("  Record length = {0:4d}, reduced to {1:4d}".format(len(pulse_list[0]), record_length))
    delta=len(pulse_list[0])-record_length
    pulse_list=pulse_list[:, int(delta/2):record_length+int(delta/2)] # pulse records are longuer than noise records (by 2)

    # Generate pulse template as average of detected pulses
    pulse_template = pulse_average(pulse_list,ignore_outlayers=True)
    if verbose:
        print("  Pulse template: ",pulse_template)

    # ############################################################
    # Comparison with xifusim data
    # ############################################################
    print("\nComparing xifusim and dre data...")
    
    # Load xifusim template
    xifusim_time,xifusim_template = np.load(file_xifusim_template)
    # Re-sample xifusim template
    xifusim_template_decimated = np.convolve(xifusim_template,1/128.*np.ones(128),mode="valid")[79::128][:1000] # Fitted best phase by hand
    
    # Determine scaling factor between xifusim and DRE units with the baseline value
    dre_template = pulse_template
    baseline_scaling = dre_template[0]/xifusim_template[0]
    PH_scaling = (dre_template[0]-dre_template.min())/(xifusim_template[0]-xifusim_template.min())
    print("  Difference between average baseline and PH scaling: {0:5.3f}%".format((PH_scaling/baseline_scaling-1)*100))
                    
    # Computing power spectra
    xifusim_PS = abs(np.fft.rfft(xifusim_template_decimated*baseline_scaling))**2
    decal=185  # Shift of DRE template to match xifusim template 
    DRE_PS = abs(np.fft.rfft(dre_template[decal:1000+decal]))**2
    PS_freq = np.fft.fftfreq(1000,6.4e-6)[:501]

    # ############################################################
    # Addition of TES noise ontop DRE noise
    # ############################################################
    print("\nComparing TES noise and DRE noise...")
    # Load TES noise spectrum
    tes_noise = fits.getdata(file_xifusim_tes_noise,"SPEC{0:04d}".format(record_length))["CSD"]*baseline_scaling
    total_noise = np.sqrt(tes_noise**2+noise_spectrum**2)

    # ############################################################
    # Compute optimal filters
    # ############################################################
    print('\nComputing optimal filter...')
    optimal_filter = compute_optimal_filter(pulse_template,noise_spectrum,energy=7.)
    #np.save('optimal_filter.npy', optimal_filter)
    print("Check optimal filter normalization: {0:5.3f}keV".format(np.dot(optimal_filter,pulse_template)))

    # Compute optimal filter, now including TES noise
    print('\nComputing optimal filter including TES noise...')
    optimal_filter_tot = compute_optimal_filter(pulse_template,total_noise,7.)
    print("Check optimal filter normalization: {0:5.3f}".format(np.dot(optimal_filter_tot,pulse_template)))

    # ############################################################
    # Do plots
    # ############################################################
    if do_plots:
        fig = plt.figure(figsize=(8, 10))    

        # Show noise spectrum
        ax1=fig.add_subplot(4,2,1)
        ax1.loglog(frequencies,noise_spectrum[1:int(record_length/2)])
        ax1.set_title('Average noise spectrum')
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel("ADU/rHz")
        for item in ([ax1.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax1.xaxis.label, ax1.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            item.set_fontsize(6)
                
        # Show pulse template
        ax2=fig.add_subplot(4,2,2)
        ax2.plot((np.arange(len(pulse_template))-200)*6.4e-3,pulse_template, label='Pulse template')
        ax2.plot((np.arange(len(pulse_template[:PREBUFF]))-200)*6.4e-3,pulse_template[:PREBUFF],'r',label='Pre-buffer')
        ax2.set_title('Pulse template')
        ax2.set_xlabel('Time [ms]')
        ax2.set_ylabel("ADU")
        ax2.set_ylim(0, 2**15)
        ax2.legend(loc='best', prop=dict(size=7))
        for item in ([ax2.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax2.xaxis.label, ax2.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(6)

        # Show difference between xifusim template and DRE template
        ax3=fig.add_subplot(4,2,3)
        ax3.plot(xifusim_template_decimated*baseline_scaling,label='xifusim')
        ax3.plot(dre_template[decal:1000+decal],label='dre')
        ax3.plot(xifusim_template_decimated*baseline_scaling - dre_template[decal:1000+decal],label='difference')
        ax3.set_title('Direct difference between both templates')
        ax3.legend(loc="best", prop=dict(size=7))
        for item in ([ax3.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax3.xaxis.label, ax3.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax3.get_xticklabels() + ax3.get_yticklabels()):
            item.set_fontsize(6)

        # Relative difference
        ax4=fig.add_subplot(4,2,4)
        ax4.plot((xifusim_template_decimated*baseline_scaling - dre_template[decal:1000+decal])/dre_template[decal:1000+decal]*100,label='difference')
        ax4.set_ylabel('Relative difference [%]')
        ax4.set_title('Relative difference between both templates')
        for item in ([ax4.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax4.xaxis.label, ax4.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax4.get_xticklabels() + ax4.get_yticklabels()):
            item.set_fontsize(6)

        # Show different power spectra
        ax5=fig.add_subplot(4,2,5)
        ax5.loglog(PS_freq,xifusim_PS,label="xifusim")
        ax5.loglog(PS_freq,DRE_PS,label="dre")
        ax5.set_xlabel("Frequency [Hz]")
        ax5.set_ylabel("PSD [AU]")
        ax5.set_title('Comparison of power spectra')
        ax5.legend(loc="best", prop=dict(size=7))
        for item in ([ax5.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax5.xaxis.label, ax5.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax5.get_xticklabels() + ax5.get_yticklabels()):
            item.set_fontsize(6)

        # Compare DAC noise with TES noise
        ax6=fig.add_subplot(4,2,6)
        ax6.loglog(frequencies,tes_noise[1:int(record_length/2)],label="TES noise")
        ax6.loglog(frequencies,noise_spectrum[1:int(record_length/2)],label="DAC noise")
        ax6.loglog(frequencies,total_noise[1:int(record_length/2)],label="Total noise")
        ax6.set_title('Average noise spectrum')
        ax6.set_title("DAC noise vs TES noise")
        ax6.legend(loc="best", prop=dict(size=7))
        for item in ([ax6.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax6.xaxis.label, ax6.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax6.get_xticklabels() + ax6.get_yticklabels()):
            item.set_fontsize(6)

        # Show optimal filter without TES noise
        ax7=fig.add_subplot(4,2,7)
        ax7.plot(optimal_filter_tot)
        ax7.set_title('Optimal filter')
        for item in ([ax7.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax7.xaxis.label, ax7.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax7.get_xticklabels() + ax7.get_yticklabels()):
            item.set_fontsize(6)

        # Show optimal filter with TES noise
        ax8=fig.add_subplot(4,2,8)
        ax8.plot(optimal_filter_tot)
        ax8.set_title('Optimal filter including TES noise')
        for item in ([ax8.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax8.xaxis.label, ax8.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax8.get_xticklabels() + ax8.get_yticklabels()):
            item.set_fontsize(6)

        fig.tight_layout()
        plt.savefig(os.path.join(plotdirname,'PLOT_E_RESOL_TEMPLATES.png'),bbox_inches='tight')

    return(optimal_filter, optimal_filter_tot)


# ############################################################
# Plot of energy resolution measurements
# ############################################################
def plot_er(NONLINEAR_FACTOR,array_to_fit1,bins1,coeffs1,axe_fit1,hist_fit1,baselines,energies,bl_correct_poly1,energies_c_bl, \
        array_to_fit2,bins2,coeffs2,axe_fit2,hist_fit2,phases,ph_correct_poly1,energies_c_ph,\
        array_to_fit3,bins3,coeffs3,axe_fit3,hist_fit3,c,tes_text,plotfilename,file_measures):
    fig = plt.figure(figsize=(8, 14))

    # defining a text to indicate the data set on the resolution plot
    file_measures_split=file_measures.split('/')
    if len(file_measures_split)==1:
        file_measures_split=file_measures.split('\\')
    txt_dirname='Data: '+file_measures_split[-1]
    
    # Show initial energy error
    ax1=fig.add_subplot(5, 1, 1)
    ax1.hist(array_to_fit1,bins=bins1,histtype='stepfilled',facecolor=c)
    ax1.plot(axe_fit1,hist_fit1,c='r',linewidth=2, label='Fit : Er={0:5.3f} x {1:5.3f} = {2:5.3f} eV' \
        .format(NONLINEAR_FACTOR, 2.355*coeffs1[2], 2.355*coeffs1[2]*NONLINEAR_FACTOR))
    ax1.legend(loc='upper left', prop=dict(size=7))
    ax1.set_title('Initial energy resolution '+tes_text)
    ax1.set_xlabel("Error [eV]")
    ax1.set_ylabel("Occurences")
    for item in ([ax1.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax1.xaxis.label, ax1.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(6)

    ax2=fig.add_subplot(5,2,3)
    ax2.plot(baselines,(energies-7)*1000,marker='.',linestyle='', c=c)
    ax2.plot(np.sort(baselines),(bl_correct_poly1(np.sort(baselines))-7)*1000,c='r',linewidth=2, label='Fit')
    ax2.legend(loc='upper left', prop=dict(size=7))
    ax2.set_title('Before baseline correction')
    ax2.set_xlabel('Baseline')
    ax2.set_ylabel('Energy - 7000 [eV]')
    for item in ([ax2.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax2.xaxis.label, ax2.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(6)
    
    #Correct for correlation
    ax3=fig.add_subplot(5,2,4)
    ax3.plot(baselines,(energies_c_bl-7)*1000,marker='.',linestyle='', c=c)
    ax3.set_title('After baseline correction')
    ax3.set_xlabel('Baseline')
    ax3.set_ylabel('Energy - 7000 [eV]')
    for item in ([ax3.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax3.xaxis.label, ax3.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(6)

    # Show energy error after baseline correction
    ax4=fig.add_subplot(5, 1, 3)
    ax4.hist(array_to_fit2,bins=bins2,histtype='stepfilled',facecolor=c)
    ax4.plot(axe_fit2,hist_fit2,c='r',linewidth=2, label='Fit : Er={0:5.3f} x {1:5.3f} = {2:5.3f} eV' \
        .format(NONLINEAR_FACTOR, 2.355*coeffs2[2], 2.355*coeffs2[2]*NONLINEAR_FACTOR))
    ax4.legend(loc='upper left', prop=dict(size=7))
    ax4.set_title('Energy resolution after baseline correction '+tes_text)
    ax4.set_xlabel("Error [eV]")
    ax4.set_ylabel("Occurences")
    for item in ([ax4.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax4.xaxis.label, ax4.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax4.get_xticklabels() + ax4.get_yticklabels()):
        item.set_fontsize(6)

    ax5=fig.add_subplot(5,2,7)
    ax5.plot(phases,(energies_c_bl-7)*1000,marker='.',linestyle='',c=c)
    ax5.plot(np.sort(phases),(ph_correct_poly1(np.sort(phases))-7)*1000,c='r',linewidth=2, label='Fit')
    ax5.legend(loc='upper left', prop=dict(size=7))
    ax5.set_title('Before phase correction')
    ax5.set_xlabel('Phase (samples)')
    ax5.set_ylabel('Energy-7000 [eV]')
    for item in ([ax5.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax5.xaxis.label, ax5.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax5.get_xticklabels() + ax5.get_yticklabels()):
        item.set_fontsize(6)
    
    ax6=fig.add_subplot(5,2,8)
    ax6.plot(phases,(energies_c_ph-7)*1000,marker='.',linestyle='', c=c)
    ax6.set_title('After phase correction')
    ax6.set_xlabel('Phase (samples)')
    ax6.set_ylabel('Energy-7000 [eV]')
    for item in ([ax6.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax6.xaxis.label, ax6.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax6.get_xticklabels() + ax6.get_yticklabels()):
        item.set_fontsize(6)

    # Show energy error after phase correction
    ax7=fig.add_subplot(5, 1, 5)
    ax7.hist(array_to_fit3,bins=bins3,histtype='stepfilled',facecolor=c, label=txt_dirname)
    ax7.plot(axe_fit3,hist_fit3,c='r',linewidth=2, label='Fit : Er={0:5.3f} x {1:5.3f} = {2:5.3f} eV' \
        .format(NONLINEAR_FACTOR, 2.355*coeffs3[2], 2.355*coeffs3[2]*NONLINEAR_FACTOR))
    ax7.legend(loc='upper left', prop=dict(size=7))
    ax7.set_title('Energy resolution after baseline and phase corrections '+tes_text)
    ax7.set_xlabel("Error [eV]")
    ax7.set_ylabel("Occurences")
    for item in ([ax7.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax7.xaxis.label, ax7.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax7.get_xticklabels() + ax7.get_yticklabels()):
        item.set_fontsize(6)

    fig.tight_layout()
    plt.savefig(plotfilename,bbox_inches='tight')


# ############################################################
# Pulse reconstruction
# ############################################################
def measure_er(file_measures, optimal_filter, optimal_filter_tot, pixeldirname, plotdirname, verbose=False, do_plots=True):
    """Perform the operations to measure the energy resolution (with and without tes noise).
    
    Arguments:
        - file_measures: fits file containing DRE pulses data
        - optimal_filter: optimal filter computed without the TES noise
        - optimal_filter_tot: optimal filter computed with the TES noise
        - pixeldirname: directory containing pixel's informations
        - plotdirname: location of plotfiles
        - verbose: if True some informations are printed (Default=False)
        - do_plots: if True a plot is done with all the intermediate results (Default=True)
        
    Returns: Nothing       
    """
    bcorr=3 # Order of polynomial baseline correction
    pcorr=8 # Order of polynomial arrival phase correction
    print('\nReconstructing pulses...')

    # Load pixel non-linearity factor from an information file
    NONLINEAR_FACTOR=get_nonlinear_factor(pixeldirname, verbose=verbose)
    print("  Loading pixel non linearity factor at 7keV: ", NONLINEAR_FACTOR)

    # Load pulse data containing the events to reconstruct
    print("  Loading measured pulse data from file ", file_measures)
    _, _, pulse_list, _=get_records_from_fits(file_measures, verbose=verbose)
    print("  Record length = {0:4d}".format(len(pulse_list[0])))
        
    # Compute first energy estimates
    energies,phases,baselines = energy_reconstruction(pulse_list, optimal_filter, PREBUFF, JITTER_MARGIN)

    coeffs1, array_to_fit1, bins1, axe_fit1, hist_fit1 = hist_and_fit((energies-7)*1000, 100)
    print("Resolution without TES noise (prior correction): {0:5.3f}eV".format(coeffs1[2]*2.355*NONLINEAR_FACTOR))

    # Apply baseline correction
    energies_c_bl, bl_correct_poly1 = apply_baseline_correction(energies, baselines, bcorr)
    coeffs2, array_to_fit2, bins2, axe_fit2, hist_fit2  = hist_and_fit((energies_c_bl-7)*1000, 100)
    print("Resolution without TES noise (after baseline correction): {0:5.3f}eV".format(coeffs2[2]*2.355*NONLINEAR_FACTOR))

    # Apply phase correction
    energies_c_ph, ph_correct_poly1 = apply_phase_correction(energies_c_bl, phases, pcorr)
    coeffs3, array_to_fit3, bins3, axe_fit3, hist_fit3  = hist_and_fit((energies_c_ph-7)*1000, 100)

    eres_notesnoise = coeffs3[2]*2.355*NONLINEAR_FACTOR
    print("Final resolution without TES noise (after phase correction): {0:5.3f}+-{1:5.3f}eV".format(eres_notesnoise,eres_notesnoise/(np.sqrt(2.*len(energies)))))

    if do_plots:
        plotfilename=os.path.join(plotdirname,'PLOT_E_RESOL_NO_TES_NOISE.png')
        plot_er(NONLINEAR_FACTOR,array_to_fit1,bins1,coeffs1,axe_fit1,hist_fit1,baselines,energies,bl_correct_poly1,energies_c_bl, \
            array_to_fit2,bins2,coeffs2,axe_fit2,hist_fit2,phases,ph_correct_poly1,energies_c_ph,\
            array_to_fit3,bins3,coeffs3,axe_fit3,hist_fit3,'g','(no TES noise)',plotfilename,file_measures)


    # ############################################################
    # Pulse reconstruction with OF also containing TES noise
    # ############################################################
    print("\nReconstruction with OF including TES noise...")
        
    # Compute first energy estimates
    energies,phases,baselines = energy_reconstruction(pulse_list, optimal_filter_tot, PREBUFF, JITTER_MARGIN)

    coeffs1, array_to_fit1, bins1, axe_fit1, hist_fit1 = hist_and_fit((energies-7)*1000, 100)
    print("Resolution without TES noise (prior correction): {0:5.3f}eV".format(coeffs1[2]*2.355*NONLINEAR_FACTOR))

    # Apply baseline correction
    energies_c_bl, bl_correct_poly1 = apply_baseline_correction(energies, baselines, bcorr)
    coeffs2, array_to_fit2, bins2, axe_fit2, hist_fit2  = hist_and_fit((energies_c_bl-7)*1000, 100)
    print("Resolution without TES noise (after baseline correction): {0:5.3f}eV".format(coeffs2[2]*2.355*NONLINEAR_FACTOR))

    # Apply phase correction
    energies_c_ph, ph_correct_poly1 = apply_phase_correction(energies_c_bl, phases, pcorr)
    coeffs3, array_to_fit3, bins3, axe_fit3, hist_fit3  = hist_and_fit((energies_c_ph-7)*1000, 100)

    eres_tesnoise = coeffs3[2]*2.355*NONLINEAR_FACTOR
    print("Final resolution with TES noise (after phase correction): {0:5.3f}+-{1:5.3f}eV".format(eres_tesnoise,eres_tesnoise/(np.sqrt(2.*len(energies)))))
        
    if do_plots:
        plotfilename=os.path.join(plotdirname,'PLOT_E_RESOL_WITH_TES_NOISE.png')
        plot_er(NONLINEAR_FACTOR,array_to_fit1,bins1,coeffs1,axe_fit1,hist_fit1,baselines,energies,bl_correct_poly1,energies_c_bl, \
            array_to_fit2,bins2,coeffs2,axe_fit2,hist_fit2,phases,ph_correct_poly1,energies_c_ph,\
            array_to_fit3,bins3,coeffs3,axe_fit3,hist_fit3,'b','(with TES noise)',plotfilename,file_measures)

    # np.save('energies.npy', energies)


# ############################################################
def ep(fulldirname, config, verbose=False):

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    general_tools.checkdir(plotdirname)
    pixeldirname = os.path.normcase("./Pixel_data_LPA75um_AR0.5/")
    file_xifusim_template = os.path.join(pixeldirname,"pulse_withBBFB.npy")
    file_xifusim_tes_noise = os.path.join(pixeldirname,"noise_spectra_bbfb_noFBDAC.fits")

    EP_filter_exist = False

    # searching data files
    f_type_deb = 15
    list_file_pulses = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[f_type_deb:]=="_mk_EP_filter_events_record.fits"]
    f_type_deb = 15
    list_file_noise = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[f_type_deb:]=="_record.fits"]
    f_type_deb = 15
    list_file_measures = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[f_type_deb:]=="_meas_E_resol_events_record.fits"]

    # Computing EP filter
    if len(list_file_pulses)==1 and len(list_file_noise)==1:
        file_pulses=os.path.join(datadirname, list_file_pulses[0])
        file_noise=os.path.join(datadirname, list_file_noise[0])
        optimal_filter, optimal_filter_tot=do_EP_filter(file_noise, file_pulses, file_xifusim_template, file_xifusim_tes_noise, plotdirname, verbose)
        EP_filter_exist=True
    else:
        print("No file available for EP processing")

    # Measuring energies
    if EP_filter_exist and len(list_file_measures)==1:
        file_measures=os.path.join(datadirname, list_file_measures[0])
        measure_er(file_measures, optimal_filter, optimal_filter_tot, pixeldirname, plotdirname, verbose)

# ############################################################

