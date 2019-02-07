def number_of_bins(array_to_bin):
    """
    Returns the optimal number of bins
    - array_to_bin (array) array needing binning
    """
    import numpy as np
    from scipy.stats import iqr
    
    return np.int(np.floor((np.max(array_to_bin)-np.min(array_to_bin))/(len(array_to_bin)**(-1/3.)*iqr(array_to_bin)*2)))

def gauss(x, a, b, c):
    """
    Creates a gaussian function with x values.
    - x is the axis (np.array)
    - a is the amplitude (np.float)
    - b is the mean value (np.float)
    - c is the standard deviation (np.float)
    """
    import numpy as np
    return a * np.exp(-(x - b)**2.0 / (2 * c**2)) 


def hist_and_fit(array_to_fit,bins,show=True, pltfilename='ER', inf=None, out=False):
    """ 
    Takes ans array of numbers and fits a gaussian on it
    - array_to_fit is the array to be fitted (np.array)
    - bins is the number of bins wished (np.int)
    - inf is an information array for the plot (list of strings)
    Returns optimal coefficients for gaussian
    Inspired from homonymous routine by P.Peille -- improved
    """
    
    import numpy as np
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size    
    mpl.rcParams['ytick.labelsize'] = label_size
    axis_font = {'fontname':'Arial', 'size':'12'} 
    
    #Remove outliers if needed
    if out:
        array_to_fit=array_to_fit[np.abs(array_to_fit-array_to_fit.mean())<5*array_to_fit.std()]
    
    hist,bins = np.histogram(array_to_fit,bins=bins) #Bin the data
    bin_centers = (bins[:-1] + bins[1:])/2 #Centers of bins
    p0=[np.max(hist)/(np.sqrt(2*np.pi)*array_to_fit.std()),array_to_fit.mean(),array_to_fit.std()] #Initial conditions
    bounds=((0, -np.inf, 0.),(np.inf, np.inf, 10*array_to_fit.std())) #/!\ sometimes may need to comment bounds and remove it from following line
    coeff, cov=curve_fit(gauss, bin_centers, hist, p0=p0, bounds=bounds, maxfev=10000000) #Do the actual fit
    
    #Plot the whole thing
    if show:
        plt.figure()
        axe_fit = np.linspace(array_to_fit.min(),array_to_fit.max(),200)
        ax = plt.gca()
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        if (inf is not None) and (np.array(inf).size!=0):
            par=inf[0]
            unit=inf[1]
            plt.xlabel(par + '(' + unit + ')', **axis_font)
        else:
            par=None
            unit='eV'
        for item in (ax.get_xticklabels()):
            item.set_rotation(45)

        hist_fit = gauss(axe_fit, *coeff)
        plt.hist(array_to_fit, bins=bins, facecolor='lightgreen', alpha=0.9, label=par)
        plt.plot(axe_fit,hist_fit,'r--',linewidth=2, label='Gaussian fit')
        plt.xlabel('Energy (eV)', **axis_font)
        plt.ylabel('Counts', **axis_font)
        plt.legend(loc='best', prop={'size':12})
        plt.title(r'{0:3d} counts, $\mu$ = {1:4.2f} eV and FWHM = {2:4.2f} eV'\
                .format(np.sum(hist), coeff[1], coeff[2]*2.35482), **axis_font)
        #plt.title(', $\mu$ = %.2f ' %(coeff[1]) + unit + ' and FWHM = %.2f ' %(coeff[2]*2.35482) + unit, **axis_font)
        #plt.show(block=False)
        plt.savefig(pltfilename+'.png', bbox_inches='tight')

    return coeff, cov    #Return fitting coefficients and covariance matrix