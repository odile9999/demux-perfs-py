import numpy as np
from numpy.fft import rfft
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import converter
from astropy.io import fits


###############################################################################
def get_pulse_shape_from_file(file):
    """Gets pulse shape from file. 
    3 file formats are considered (fits, npy and txt)
    
    Arguments:
        - filename
        
    Returns: pulse_shape
        - t: numpy array with time of each sample
        - pulse: numpy array with value of each sample
    """

    file_extention=file.split('.')[-1]
    if file_extention=='txt':
        f = open(file, "r")
        t=[]
        pulse=[]
        for x in f:
            t.append(float(x.split(' ')[0])) 
            pulse.append(float(x.split(' ')[1]))
        t=np.array(t)
        pulse=np.array(pulse)
    elif file_extention=='npy':
        tab=np.load(file)
        t=tab[0]
        pulse=tab[1].max()-tab[1]   ## Put the data in the correct orientation
    elif file_extention=='fits':
        hdul = fits.open(file)
        template=hdul['TEMPLATE'].data
        t=template['TIME']
        pulse=template['CURRENT'].max()-template['CURRENT']
    else :
        print("Error, wrong file type")
    return(t-t[0], pulse)

###############################################################################
def mk_LUT_code(pulse_file, ratio_nbits, output_file, fs_firmware=19.53e6):
    """Gets pulse shape from file, does the sampling to fit in the firmware
    pulse geerator and stores the resulting data in a text file.
    
    Arguments:
        - filename: name of the input data
        - ratio_nbits: scaling factor to be applied on the time array by the
                       firmware. This factor is needed to use optimally the 
                       fake pulse memory.
        - output_file: name of the output file 
        - fs_firmware: computing frequency of the firmware (default=19.53e6)
        
    Returns: pulse_shape
        - t: numpy array with time of each sample (after resampling)
        - pulse: numpy array with value of each sample (after resampling)
    """
    # reading pulse and time data from file
    t, pulse=get_pulse_shape_from_file(pulse_file)

    # re-sampling pulse data
    func = interp1d(t, pulse, kind='cubic')
    fs_model=1./(t[1]-t[0])
    print("Fs_model={0:8.2f} Hz, Fs_firmware={1:8.2f} Hz".format(fs_model, fs_firmware))
    ratio=2**ratio_nbits
    f=fs_firmware/ratio
    rom_size_nbits=11
    rom_size=2**rom_size_nbits
    t2=np.arange(rom_size)/f
    table=func(t2)

    # scaling pulse amplitude
    max_val=2**16-1
    table_pulse=(table*max_val/table.max()).astype(int)

    # computing delta table for interpolation
    table_diff=np.append(table_pulse[1:]-table_pulse[:-1], 0)

    # convertion to 2's complement values
    pulse_nbits = int(np.ceil(np.log10(abs(table_pulse).max())/np.log10(2)))
    if pulse_nbits<=16:
        pulse_nbits=16
    else:
        print("problem number of bits in the pulse is too high!")
    delta_nbits = int(np.ceil(np.log10(2*abs(table_diff).max())/np.log10(2)))
    if delta_nbits<=16:
        delta_nbits=16
    else:
        print("problem number of bits in the delta is too high!")

    f=open(output_file,"w")
    timescale=ratio_nbits-4
    f.write("{0:2d}\n".format(timescale))
    for i in range(rom_size):
        pulse_part=converter.switch_bin2hexa(converter.dec2natbin(table_pulse[i], pulse_nbits))
        delta_part=converter.switch_bin2hexa(converter.dec2cad(table_diff[i], delta_nbits))
        f.write('0x'+pulse_part+delta_part+'\n')
    f.close()

    return(t2, table_pulse)

###############################################################################

dir_in=os.path.normcase('./files_in/')
dir_out=os.path.normcase('./files_out/')
CBE_f_in_txt=os.path.join(dir_in,'pulse_75um_AR0_5_withoutBBFB.npy')
CBE_f_in_fits=os.path.join(dir_in,'XIFU_INST_XIFUSIM_20191205-164022_PP_PARAM_CONFXIFUSIM.fits')
CBE_f_out=os.path.join(dir_out,'./fake_pulse_CBE.txt')
banc_f_in=os.path.join(dir_in,'./pulse_banc_test.txt')
banc_f_out=os.path.join(dir_out,'./fake_pulse_banc.txt')

t1, pulse1=mk_LUT_code(CBE_f_in_txt, 6, CBE_f_out)
t2, pulse2=mk_LUT_code(CBE_f_in_fits, 6, CBE_f_out)
t3, pulse3=mk_LUT_code(banc_f_in, 5, banc_f_out)

plt.plot(t1-t2)
plt.show()