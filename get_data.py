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
def readfile(dumpfilename, quiet=True):
    r"""
        This function reads data from a DRE dump file, and returns 2 arrays
        (format <h).
        The file header (first 32-bit) word is removed from the data

        Parameters
        ----------
        dumpfilename : string
        The name of the dump file (with the path and the extension)

        quiet : boolean
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
    #header21=header2-header24*2**12-header23*2**8-header22*2**4 not used
    dumptype=header22
    if not quiet:
        print('  Dump type is: ' + dumpstr(dumptype))
    
    data = np.resize(data, (len(data)//2, 2))

    return(data, dumptype)

# -----------------------------------------------------------------------------
def read_iq(filename):
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
    
    dada_ch0 = tab_i[:, 0]
    dada_ch1 = tab_i[:, npix+2]
    
    chan0_id = tab_q[:, 0]
    chan1_id = tab_q[:, npix+2]
    
    chan0_i = tab_i[:, 2:2+npix]
    chan1_i = tab_i[:, 2+npix+2:]
    
    chan0_q = tab_q[:, 2:2+npix]
    chan1_q = tab_q[:, 2+npix+2:]

    flag_error = check_data(dada_ch0, dada_ch1, chan0_id, chan1_id)

    return(chan0_i, chan0_q, chan1_i, chan1_q, flag_error)

# -----------------------------------------------------------------------------
def check_data(dada_ch0, dada_ch1, chan0_id, chan1_id):
    r"""
        This function reads IQ data from a DRE IQ file.

        Parameters
        ----------
        dada_ch0, dada_ch1, chan0_id, chan1_id : array type
        dada flags and channel ids for channels 0 and 1

        Returns
        -------
        FLAG_ERROR : boolean
        True if the data format is not correct
                    
        """
    flag_error = False
    dada = -9510
    ch0_id = 10880
    ch1_id = 10881
    check_dada0 = np.where(dada_ch0 != dada)
    check_dada1 = np.where(dada_ch1 != dada)
    check_ch0 = np.where(chan0_id != ch0_id)
    check_ch1 = np.where(chan1_id != ch1_id)

    if len(check_dada0[0]) > 0 or len(check_dada1[0]) > 0 \
        or len(check_ch0[0]) > 0 or len(check_ch1[0]) > 0:
        print(" Error! Problem in the data set !!!!")
        flag_error = True        

    return(flag_error)

# -----------------------------------------------------------------------------
def read_events(eventsfilename, backup_version):
    r"""
        This function reads events binary files

        Parameters
        ----------
        eventsfilename : string
        The name of the dump file (with the path and the extension)

        backup_version : number
        reference of the backup_version version

        Returns
        -------
        timestamp : array
        photon start time (s)

        channelId : array
        channel which measured the photon

        pixelId : array
        pixel which measured the photon

        energy : array
        energy of the photon (eV)

        offset : array
        value of the baseline measured before the event
        """
    
    dt1=np.dtype([('timestamp', np.float), \
                 ('pixelId', np.int8), \
                 ('energy', np.float32), \
                 ('offset', np.float32)])

    dt2=np.dtype([('timestamp', np.float), \
                 ('channelId', np.int8), \
                 ('pixelId', np.int8), \
                 ('energy', np.float32), \
                 ('offset', np.float32)])
                                  
    if backup_version < 1:
        fdat=open(eventsfilename, 'rb')
        event_list=np.fromfile(fdat, dtype=dt1)
        channel_id = np.zeros(len(event_list[:]['timestamp']))
        fdat.close()
        return(event_list[:]['timestamp'], channel_id, event_list[:]['pixelId'], event_list[:]['energy'], event_list[:]['offset'])
    else:
        fdat=open(eventsfilename, 'rb')
        event_list=np.fromfile(fdat, dtype=dt2)
        fdat.close()
        return(event_list[:]['timestamp'], event_list[:]['channelId'], event_list[:]['pixelId'], event_list[:]['energy'], event_list[:]['offset'])
    
# -----------------------------------------------------------------------------
