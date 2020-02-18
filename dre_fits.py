import params
import numpy as np
from astropy.io import fits


# -----------------------------------------------------------------------------
def read_fits(filename, verbose=False):
    r"""
        This function reads data from a DRE dump file (fits format)
        It returns 2 arrays.
        The file header (first 32-bit) word is removed from the data

        Parameters
        ----------
        filename : string
        The name of the dump file (with the path and the extension)

        verbose : boolean
        If True informations are written by the routine (default=False)

        Returns
        -------
        ftype : string
        describes the file type

        header : array type
        header of the fits file.

        data : array type
        data of the fits file.
            
        """
    print("Reading data from fits file...")
    hdul=fits.open(filename)
    header=hdul[1].header
    data=hdul[1].data
    ftype=header['FILETYPE']

    if verbose:
        print("  Informations of FITS file:")
        print("    Date:      ", header['DATE'])
        print("    Origin:    ", header['ORIGIN'])
        print("    Project:   ", header['INSTRUME'])
        print("    File type: ", ftype)

    return(ftype, header, data)

# -----------------------------------------------------------------------------
def read_iq(filename, verbose=False):
    r"""
        This function reads IQ data from a DRE fits file.
        (this corresponds to standard observation data)

        Parameters
        ----------
        filename : string
        The name of the dump file (with the path and the extension)

        verbose : boolean
        If True informations are written by the routine (default=False)

        Returns
        -------
        i, q : array type of floats
        contains the values of the I and Q science data.
                    
        """
    ftype, header, data = read_fits(filename, verbose)

    if ftype != 'IQ':
        raise ValueError('Wrong file type!')

    nval=len(data['timestamp'][0])
    
    i=np.zeros((nval, params.npix))
    q=np.zeros((nval, params.npix))
    for pix in range(params.npix):
        if header['PIXEL_{0:d}'.format(pix)]==1:  # pixel is active
            i[:,pix]=data['I{0:d}'.format(pix)][0]
            q[:,pix]=data['Q{0:d}'.format(pix)][0]

    return(i.astype('float'), q.astype('float'))

# -----------------------------------------------------------------------------
def read_dump(filename, checktype='NO-CHECK', verbose=False):
    r"""
        This function reads dumps data from a DRE file.
        (this corresponds to diagnostic data)

        Parameters
        ----------
        filename: string
        The name of the dump file (with the path and the extension)

        checktype: string
        The expected file type. This file type is compared with the file 
        type wich is indicated in the fits header. Default is NO-CHECK.

        verbose: boolean
        If True informations are written by the routine (default=False)

        Returns
        -------
        field1, field2: array type of floats
        contains the values of the two data streams and the type of dump.

        ftype: string                
        dump type readout from the header of the fits file.
        """
    ftype, header, data = read_fits(filename, verbose)

    if checktype not in [ftype, 'NO-CHECK']:
        raise ValueError('Wrong file type!')

    return(data['Field1'][0].astype('float'), data['Field2'][0].astype('float'), ftype)

# -----------------------------------------------------------------------------
def read_records(filename, verbose=False):
    r"""
        This function reads sample records from a fits file.

        Parameters:
        ----------
        filename: string
        The name of the fits file (with the path and the extension)

        verbose: boolean
        If True informations are written by the routine (default=False)

        Returns
        -------
        Returns: arrays containing the channel id, the pixel id, the records of module, the time stamp
    """

    hdul = fits.open(filename)
    data=hdul[1].data
    t=data['Timestamp']
    chid=data['channelNum']
    pixid=data['pixelNum']
    i=data['i'].astype("float")
    q=data['q'].astype("float")
    module=np.sqrt(i**2+q**2)
    if verbose:
        print("  Informations of FITS file:")
        print("    Date:    ", hdul[1].header['DATE'])
        print("    Origin:  ", hdul[1].header['ORIGIN'])
        print("    Project: ", hdul[1].header['INSTRUME'])
        print("    Number of records: ", len(i))
        print("    Length of records: ", int(len(i[0])))
    
    return(chid, pixid, module, t)

# -----------------------------------------------------------------------------



# filename='test_IQ.fits'
# i,q=read_iq(filename, True)

# filename2='20190411_101402_0000_IN-BIA_PULSE.fits'
# f1,f2,dumptype=read_dump(filename2, True)

# filename3='20190411_101305_0011_IQ-TST_Science-Data.fits'
# itst,qtst=read_iqtst(filename3, True)
# modtst=np.sqrt(itst.astype(float)**2+qtst.astype(float)**2)