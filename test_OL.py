import numpy as np
from numpy.fft import rfft
import os
import matplotlib.pyplot as plt
import get_data, general_tools


# Function to apply a low pass filter
def smooth(x,window_len=11,window='hanning'):

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


def test_OL(fulldirname, config):

    datadirname = os.path.join(fulldirname, config['dir_data'])
    L=2**18

    # Unfiltered data
    fich='20190802_155718_0040_IQ-TST_Open_Loop.dat'
    data, _ = get_data.readfile(os.path.join(datadirname, fich))
    i=data[1:L,0]
    q=data[1:L,1]
    modulus = np.sqrt(i.astype('float')**2 + q.astype('float')**2)

    # Filtered data
    fich='20190802_155720_0040_IQ-ALL_Open_Loop.dat'
    i_filt, q_filt, _, _, _ = get_data.read_iq(os.path.join(datadirname, fich))
    i_filt=i_filt[1:L,0]
    q_filt=q_filt[1:L,0]
    modulus_filt = np.sqrt(i_filt.astype('float')**2 + q_filt.astype('float')**2)


    i_filt=smooth(i,40)
    q_filt=smooth(q,40)
    modulus_filt = np.sqrt(i_filt**2 + q_filt**2)


    print("I : {0:6d} values".format(len(i)))
    print("Q : {0:6d} values".format(len(q)))

    fig=plt.figure(1,(8,6))
    ax1=fig.add_subplot(3,1,1)
    ax1.plot(i, label='i')
    ax1.plot(q, label='q')
    ax1.plot(modulus, label='module')
    ax1.legend(loc='best')

    k=128
    ax2=fig.add_subplot(3,1,2)
    ax2.plot(i[:k], label='i')
    ax2.plot(q[:k], label='q')
    ax2.plot(modulus[:k], label='module')
    ax2.legend(loc='best')


    ax3=fig.add_subplot(3,1,3)
    ax3.plot(i_filt[80:L-80], label='i-filtered')
    ax3.plot(q_filt[80:L-80], label='q-filtered')
    ax3.plot(modulus_filt[80:L-80], label='module-filtered')
    ax3.legend(loc='best')

    fig.tight_layout()

    plt.show()