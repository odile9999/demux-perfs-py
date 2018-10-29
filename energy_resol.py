# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:28:01 2018

@author: Paul
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from numpy.linalg import inv
from DEMUX_TBox.get_data import readIQ
from scipy.optimize import curve_fit


#------------------------------------------------------------------------------
def get_test_pix(filename):

    pix_test=40
    Chan0_i, Chan0_q, Chan1_i, Chan1_q, FLAG_ERROR = readIQ(filename)
    modulus = \
        np.sqrt(Chan0_i[:,pix_test].astype('float')**2 + Chan0_q[:,pix_test].astype('float')**2)

    "Sampling frequency"
    fs=20e6/2**7
    t=np.arange(len(modulus))/fs

    return(t, modulus)

#------------------------------------------------------------------------------
def meas_energy_r(fulldirname, config):

    N=2048
    seuil=70
    ind=3
    fs=20e6*2**7

    datadirname = os.path.join(fulldirname, config['dir_data'])
    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)

    pltfilename = os.path.join(plotdirname, "PLOT_ENERGY-RESOL")

    # Reading calibration data from files
    i_type_deb, i_type_fin = 27, 36
    
    noise_name = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat'\
                and f[i_type_deb:i_type_fin]=="_ER-Noise"]
    t, noise = get_test_pix(os.path.join(datadirname, noise_name[0]))

    pulse_name = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat'\
                and f[i_type_deb:i_type_fin]=="_ER-Pulse"]
    t, pulse = get_test_pix(os.path.join(datadirname, pulse_name[0]))


    offset=np.mean(noise)
    noise=offset-noise
    pulse=offset-pulse

    pulses=[]
    pulse_spectrum=np.zeros(N)
    pulse_phase=np.zeros(N)
    Nnoise=noise.__len__()//N
    noise_spectrum=np.zeros(N)
    RI=np.zeros(N)

    # Pulses detection and pulse spectrum computation      
    j=1
    i=0
    k=0
    freq=np.arange(N)*fs/N
    fig = plt.figure(figsize=(9, 7))
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_xlabel(r'Time (samples)')
    ax1.set_ylabel(r'Amplitude (A.U.)')
    ax1.set_title(r'Calibration pulses')

    while i+j+N<pulse.__len__():
        if j==0 and pulse[i]>seuil:
            l=1
            while pulse[i+l]>seuil:
                l=l+1
            i=i+l
        else:
            if pulse[i+j]<seuil:
                j=j+1
            else:
                pulses=np.append(pulses,pulse[i+j-ind:i+j-ind+N],axis=0)
                pulse_spectrum=pulse_spectrum+abs(np.fft.fft(np.array(pulse[i+j-ind:i+j-ind+N])))
                ax1.plot(pulse[i+j-ind:i+j-ind+N])
                if k==5:
                    pulse_phase=np.angle(np.fft.fft(pulse[i+j-ind:i+j-ind+N]))
                i=i+N+j-ind
                j=0
                k=k+1
    n_pulse=k-1
    pulses=np.resize(pulses,(n_pulse,N))

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.semilogx(freq, pulse_spectrum)
    ax2.set_xlabel(r'Frequency (Hz)')
    ax2.set_ylabel(r'Amplitude (A.U.)')
    ax2.set_title(r'Mean pulses spectrum')
    ax1.text(0.5, 0.9, "{0:4d} pulses detected".format(n_pulse), color='b', transform=ax1.transAxes)

    fig.tight_layout()

    # Noise spectrum computation      
    for i in range(0,Nnoise):
        noise_spectrum=noise_spectrum+abs(np.fft.fft(noise[i*N:(i+1)*N]))

    # RI : Réponse impulsionnelle
    RI=np.real(np.fft.ifft(pulse_spectrum/noise_spectrum*np.exp(1j*pulse_phase)))

    ax = fig.add_subplot(2, 2, 4)
    ax.semilogx(freq,20*np.log10(noise_spectrum))
    ax.set_xlabel(r'Frequency (Hz)')
    ax.set_ylabel(r'Amplitude (A.U.)')
    ax.set_title(r'Mean noise spectrum')
    fig.tight_layout()

    ax = fig.add_subplot(2, 2, 3)
    ax.plot(RI)
    ax.set_xlabel(r'Time (samples)')
    ax.set_ylabel(r'Amplitude (A.U.)')
    ax.set_title(r'Impulse response')
    fig.tight_layout()
    plt.savefig(pltfilename+'_1.png', bbox_inches='tight')

    f=0
    for i in range(0,n_pulse):
        f=f+np.dot(RI,pulses[i,:])

    # scaling factor : f
    f=f/n_pulse

    data=pulse

    Np=3
    X=np.zeros(Np)
    X=np.resize(X,(Np,3))
    for i in range(0,Np):
        for j in range(0,3):
            X[i][j]=(i-Np//2)**(2-j)
            
    Z=inv(np.dot(np.transpose(X),X))
    V=np.zeros(Np)
    E=np.zeros(n_pulse)
    t0=np.zeros(n_pulse)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(1, 1, 1)
    j=1
    i=0
    k=0
    while i+j+N<data.__len__() and k<n_pulse:
        if j==0 and data[i]>seuil: # In case of too close pulses 
            l=1
            while data[i+l]>seuil:
                l=l+1
            i=i+l
        else:
            if data[i+j]<seuil:
                j=j+1
            else:
                for m in range(0,Np):  # Np (=3) considered timings (photon wrt sampling)
                    V[m]=np.dot(RI,data[i+j-ind+(-Np//2+m):i+j-ind+N+(-Np//2+m)])
                ax.plot(V)
                P=np.dot(np.dot(Z,np.transpose(X)),V)
                # t(Emax) and Emax for an order 2 fit
                t0[k]=-P[1]/(2*P[0])
                E[k]=(P[2]-P[1]**2/(2*P[0]))*(7000/f)
                i=i+N+j-ind
                j=0
                k=k+1
    fig.tight_layout()
    plt.savefig(pltfilename+'_3.png', bbox_inches='tight')


    ########################                       
    # "Banana fit"
    T=[]
    T=np.resize(T,(t0.__len__(),3))
    T[:,0]=t0**2
    T[:,1]=t0
    T[:,2]=np.ones(t0.__len__())
    pfit=np.dot(np.dot(inv(np.dot(np.transpose(T),T)),np.transpose(T)),np.transpose(E))


    ##################################
    # Processing des données observées
    pulse = noise = 0  # to free some memory space
    # Reading "observation" data from files    
    observ_name = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]=='.dat'\
                and f[i_type_deb:i_type_fin]=="_ER-Obser"]
    t, observ = get_test_pix(os.path.join(datadirname, observ_name[0]))
    observ=offset-observ

    data = observ
    # Counting the number of events
    j=1
    i=0
    k=0
    while i+j+N<data.__len__():
        if j==0 and data[i]>seuil:
            l=1
            while data[i+l]>seuil:
                l=l+1
            i=i+l
        else:
            if data[i+j]<seuil:
                j=j+1
            else:
                i=i+N+j-ind
                j=0
                k=k+1
    n_pulse = k-1

    # Processing the events
    Ecor1=np.zeros(n_pulse)
    Ecor2=np.zeros(n_pulse)
    t0=np.zeros(n_pulse)
    j=1
    i=0
    k=0
    while i+j+N<data.__len__() and k<n_pulse:
        if j==0 and data[i]>seuil:
            l=1
            while data[i+l]>seuil:
                l=l+1
            i=i+l
        else:
            if data[i+j]<seuil:
                j=j+1
            else:
                for m in range(0,Np):
                    V[m]=np.dot(RI,data[i+j-ind+(-Np//2+m):i+j-ind+N+(-Np//2+m)])
                P=np.dot(np.dot(Z,np.transpose(X)),V)
                t0[k]=-P[1]/(2*P[0])

                # Banana correction
                fph=(pfit[0]*t0[k]**2+pfit[1]*t0[k]+pfit[2])/7000
                Ecor1[k]=(P[2]-P[1]**2/(2*P[0]))*(7000/f)
                Ecor2[k]=Ecor1[k]/fph
                
                i=i+N+j-ind
                j=0
                k=k+1

    fig = plt.figure(figsize=(9, 7))
    nbins = 30
    tmin, tmax = -1, 3
    Emin = np.floor(np.min(np.concatenate((Ecor1, Ecor2))))
    Emax = np.floor(np.max(np.concatenate((Ecor1, Ecor2))))+1
    Emaxrange = n_pulse / 8

    ax3 = fig.add_subplot(2, 2, 1)
    ax3.plot(t0, Ecor1, ".", color='g', alpha=0.75, label='Simple sampling phase correction')
    ax3.axis([tmin, tmax, Emin, Emax])
    ax3.set_xlabel(r'Time (samples)')
    ax3.set_ylabel(r'Energy (Ev)')
    ax3.legend(loc='best')

    ax4 = fig.add_subplot(2, 2, 2)
    #ax4.hist(Ecor1, nbins)
    n, bins, patches = ax4.hist(Ecor1, nbins, density=False, facecolor='g', alpha=0.75)
    sum1 = np.sum(n)/nbins
    max_gauss1 = gauss_fit(ax4, Ecor1, sum1, n_pulse)
    ax4.axis([Emin, Emax, 0, Emaxrange])
    ax4.set_ylabel(r'Counts')
    ax4.set_xlabel(r'Energy (eV)')

    ax5 = fig.add_subplot(2, 2, 3)
    ax5.plot(t0, Ecor2, ".", color='b', alpha=0.75, label='Full sampling phase correction')
    ax5.axis([tmin, tmax, Emin, Emax])
    ax5.set_xlabel(r'Time (samples)')
    ax5.set_ylabel(r'Energy (Ev)')
    ax5.legend(loc='best')

    ax6 = fig.add_subplot(2, 2, 4)
    #ax6.hist(Ecor2, nbins)
    n, bins, patches = ax6.hist(Ecor2, nbins, density=False, facecolor='b', alpha=0.75)
    sum1 = np.sum(n)/nbins
    max_gauss2 = gauss_fit(ax6, Ecor2, sum1, n_pulse)
    ax6.set_ylabel(r'Counts')
    ax6.set_xlabel(r'Energy (eV)')
    ax6.axis([Emin, Emax, 0, Emaxrange])

    max_gauss = max(max_gauss1, max_gauss2)
    ax4.set_ylim = ([0, 1.2*max_gauss])
    ax6.set_ylim = ([0, 1.2*max_gauss])

    fig.tight_layout()
    plt.savefig(pltfilename+'_2.png', bbox_inches='tight')


# -----------------------------------------------------------------------------
# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# -----------------------------------------------------------------------------
def gauss_fit(ax, E, sum1, counts):
    hist, bin_edges = np.histogram(E, density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

    min, max, npts = np.min(E), np.max(E), 100
    xE = (np.arange(npts)-npts/2) * (max-min)/npts + (max+min)/2
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [10., 7000., 1.]
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    # Get the fitted curve
    hist_fit = gauss(xE, *coeff)
    sum2 = np.sum(hist_fit)/npts
    f=sum1/sum2
    fwhm = 2.355*coeff[2]
    ax.plot(xE, hist_fit*f, '--r', label='Fitted data', linewidth=2)
    ax.text(0.05, 0.9, 'FWHM={0:4.2f}eV'.format(fwhm), transform=ax.transAxes)
    ax.text(0.05, 0.85, 'Nb counts={0:4d}'.format(counts), transform=ax.transAxes)
    return(np.max(hist_fit*f))

# -----------------------------------------------------------------------------
