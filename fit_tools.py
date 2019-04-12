def number_of_bins(array_to_bin):
    """
    Returns the optimal number of bins
    - array_to_bin (array) array needing binning
    """
    import numpy as np
    from scipy.stats import iqr    
    return np.int(np.floor((array_to_bin.max()-array_to_bin.min())/(len(array_to_bin)**(-1/3.)*iqr(array_to_bin)*2)))

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

def l_gauss(x, e, a):
    """
    Returns the line function 
    - x is an array containing FW, coef, cen, base
    - e is an array with the energy
    - a is an array with the data to fit
    - name is the line name to import right data (string)
    - bool_cont (boolean) if noise is done or not
    Returns function for lsq or other minimization algorithms
    """   
    import numpy as np

    weights=np.ones(len(a))
    weights[a>0]=1/np.sqrt(a[a>0])
    return np.sum(((gauss(e, x[0], x[1], x[2])-a)/weights)**2)

def gauss_fit(array_to_fit, bins, show=True, pltfilename='ER', inf=None):
    """ 
    Takes ans array of numbers and fits a gaussian on it
    - array_to_fit is the array to be fitted (np.array)
    - bins is the number of bins wished (np.int)
    - inf is an information array for the plot (list of strings)
    Returns optimal coefficients for gaussian
    """
    
    import numpy as np
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    nl_at_7kev = 1.361 # NL of LPA75um pixels at 7keV
    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size    
    mpl.rcParams['ytick.labelsize'] = label_size
    axis_font = {'fontname':'Arial', 'size':'12'} 
    
    hist,bins = np.histogram(array_to_fit,bins=bins)
    bin_centers = (bins[:-1] + bins[1:])/2
    p0=[1.001*np.max(hist), np.average(array_to_fit), np.std(array_to_fit)]
    res = minimize(l_gauss, p0, args=(bin_centers, hist), method='Powell', 
                        tol=1e-20, options={'maxiter':1000000, 'ftol': 1e-20})

    coeff=res.x

    if res.success==False:
        print("Fit dit not converge")
    
    #Plot the whole thing
    if show:
        plt.figure(figsize=(8, 7))
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
        plt.title(r'{0:3d} counts, $\mu$ = {1:4.2f} eV and FWHM = {2:4.3f}x{3:4.2f} = {4:4.2f} eV'\
                .format(np.sum(hist), coeff[1], nl_at_7kev, abs(coeff[2])*2.35482, nl_at_7kev*abs(coeff[2])*2.35482), **axis_font)
        plt.savefig(pltfilename+'.png', bbox_inches='tight')

def low_pass(f, k, pi2rc):
    """
    Creates a low pass filter transfer function with f values.
    - f is the axis (np.array)
    - k is the DC level (np.float)
    - pi2rc is equal to the value of 2xPixRC for the equivalent RC low pass filter (np.float)
    """
    import numpy as np
    return k/np.sqrt(1.+(f*pi2rc)**2)

def l_low_pass(x, f, a):
    """
    Returns the line function 
    - x is an array containing the function parameters
    - f is an array with the frequencies
    - a is an array with the data to fit
    Returns function for lsq or other minimization algorithms
    """   
    import numpy as np

    weights=np.ones(len(a))
    weights[a>0]=1/np.sqrt(a[a>0])
    return np.sum(((low_pass(f, x[0], x[1])-a)/weights)**2)

def low_pass_fit(freqs, array_to_fit):
    """ 
    Takes ans array of numbers and fits a low pass filter function on it
    - array_to_fit is the array to be fitted (np.array)
    Returns optimal coefficients for the fit
    """
    import numpy as np
    from scipy.optimize import minimize

    p0=[array_to_fit[0], 1/(15e3)]
    res = minimize(l_low_pass, p0, args=(freqs, array_to_fit), method='Powell', 
                        tol=1e-20, options={'maxiter':1000000, 'ftol': 1e-20})

    if res.success==False:
        print("Fit dit not converge")

    return res.x
