import numpy as np

# -----------------------------------------------------------------------------
def dumpstr(dumptype):
    r"""
        This function translates the dumptype into a string

        Parameters
        ----------
        dumptype : number
        The type of dump

        Returns
        -------
        dumptype_str : string
        the dumptype as a string.
        
        """
    # Checking the dump type
    if dumptype==0:
        dumptype_str = "Channel 0, input and bias"
    elif dumptype==4:
        dumptype_str = "Channel 1, input and bias"
    elif dumptype==1:
        dumptype_str = "Channel 0, input and feedback"
    elif dumptype==5:
        dumptype_str = "Channel 1, input and feedback"        
    elif dumptype==2:
        dumptype_str = "Channel 0, bias and feedback"        
    elif dumptype==6:
        dumptype_str = "Channel 1, bias and feedback"        
    elif dumptype==8:
        dumptype_str = "IQ data"                
    elif dumptype==9 or dumptype==10:
        dumptype_str = "test pixel IQ data"                
    elif dumptype==15:
        dumptype_str = "32-bit counter"                
    else:
        raise ValueError('Wrong dump type')
    
    return(dumptype_str)


# -----------------------------------------------------------------------------
def readfile(dumpfilename, Quiet=True):
    r"""
        This function reads data from a DRE dump file, and returns 2 arrays
        (format <h).
        The file header (first 32-bit) word is removed from the data

        Parameters
        ----------
        dumpfilename : string
        The name of the dump file (with the path and the extension)

        Quiet : boolean
        Defines if text info shall be written by the routine

        Returns
        -------
        data : array type
        contains the values of the file (format int16).
            
        dumptype : number 'uint16'
        Identifier for the type of dump.
        
        """
    
    fdat=open(dumpfilename, 'rb')
    data=np.fromfile(fdat, dtype='<h')
    fdat.close()
    
    DADA=-9510      # 0xDADA interpreted as int16
    if data[0] != DADA:
        raise ValueError('Problem with file format!')
    header2=data[1].astype('uint16')
    header24=int(header2/2**12)
    header23=int((header2-header24*2**12)/2**8)
    header22=int((header2-header24*2**12-header23*2**8)/2**4)
    #header21=header2-header24*2**12-header23*2**8-header22*2**4
    dumptype=header22
    if not Quiet:
        print('  Dump type is: ' + dumpstr(dumptype))
    
    data = np.resize(data, (len(data)//2, 2))

    return(data, dumptype)


# -----------------------------------------------------------------------------
def readIQ(filename):
    r"""
        This function reads IQ data from a DRE IQ file.
        (this corresponds to standard observation data)

        Parameters
        ----------
        filename : string
        The name of the dump file (with the path and the extension)

        Returns
        -------
        Chan0_i, Chan0_q, Chan1_i, Chan1_q : array type
        contains the values of the I and Q science data (format int16).

        FLAG_ERROR : boolean
        if FLAG_ERROR is True the file format is incorrect.
                    
        """

    data, dumptype = readfile(filename)
    if dumptype != 8:
        raise ValueError('Wrong dumptype')

    # decommutation des donnees
    npix = 41
    
    npts = len(data[:,0]) // (2*(npix+2))        
    tab_i = np.resize(data[:,0], (npts, 2*(npix+2)))
    tab_q = np.resize(data[:,1], (npts, 2*(npix+2)))
    
    DADA_ch0 = tab_i[:, 0]
    DADA_ch1 = tab_i[:, npix+2]
    
    Chan0_ID = tab_q[:, 0]
    Chan1_ID = tab_q[:, npix+2]
    
    Chan0_i = tab_i[:, 2:2+npix]
    Chan1_i = tab_i[:, 2+npix+2:]
    
    Chan0_q = tab_q[:, 2:2+npix]
    Chan1_q = tab_q[:, 2+npix+2:]

    FLAG_ERROR = check_data(DADA_ch0, DADA_ch1, Chan0_ID, Chan1_ID)

    return(Chan0_i, Chan0_q, Chan1_i, Chan1_q, FLAG_ERROR)


# -----------------------------------------------------------------------------
def check_data(DADA_ch0, DADA_ch1, Chan0_ID, Chan1_ID):
    r"""
        This function reads IQ data from a DRE IQ file.

        Parameters
        ----------
        DADA_ch0, DADA_ch1, Chan0_ID, Chan1_ID : array type
        DADA flags and channel ids for channels 0 and 1

        Returns
        -------
        FLAG_ERROR : boolean
        True if the data format is not correct
                    
        """
    FLAG_ERROR = False
    DADA = -9510
    CH0_ID = 10880
    CH1_ID = 10881
    check_DADA0 = np.where(DADA_ch0 != DADA)
    check_DADA1 = np.where(DADA_ch1 != DADA)
    check_CH0 = np.where(Chan0_ID != CH0_ID)
    check_CH1 = np.where(Chan1_ID != CH1_ID)

    if len(check_DADA0[0]) > 0 or len(check_DADA1[0]) > 0 \
        or len(check_CH0[0]) > 0 or len(check_CH1[0]) > 0:
        print(" Error! Problem in the data set !!!!")
        FLAG_ERROR = True        

    return(FLAG_ERROR)

# -----------------------------------------------------------------------------
