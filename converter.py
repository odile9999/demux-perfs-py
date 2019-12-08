
###############################################################################
def switch_bin2hexa_digit(digit):
    switcher={
    "0000": "0",
    "0001": "1",
    "0010": "2",  
    "0011": "3", 
    "0100": "4",   
    "0101": "5",   
    "0110": "6",   
    "0111": "7",   
    "1000": "8",   
    "1001": "9",   
    "1010": "A",   
    "1011": "B",   
    "1100": "C",   
    "1101": "D",   
    "1110": "E",   
    "1111": "F"  
    }
    return(switcher[digit])

###############################################################################
def switch_bin2hexa(bin):
    r"""
    This function makes the conversion from binary ascii to hexadecimal ascci.
    """

    nbits=len(bin)
    offset=0
    hexa=''
    while offset<nbits:
        hexa=hexa+switch_bin2hexa_digit(bin[offset:offset+4])
        offset+=4
    return(hexa)

###############################################################################
def dec2cad(dec, n):
    r"""
    This function computes the 2s complement binary value of an integer.
    
    Parameters
    ----------
    dec : number
        The decimal input integer to be converted
    n : number
        The number of bits of the output string

    Output
    ------
    bin_str : string
        The binary value in a string format

    Raises
    ------
    ValueError
        When n is too small to do the conversion

    See Also
    --------
    num : string to number conversion

    Examples
    --------
    >>> dec2cad(7, 8)
    '00000111'

    >>> dec2cad(-7, 8)
    '11111001'

    >>> dec2cad(7, 2)
    ValueError: Requested size is too small for this value

    """
    bin_str = ''
    if dec >= 0:
        while n > 0:
            bin_str = str(dec % 2) + bin_str
            dec >>= 1
            n += -1
    else:
        dec = abs(dec + 1)
        while n > 0:
            rev_str = ('1', '0')
            bin_str = rev_str[dec % 2] + bin_str
            dec >>= 1
            n += -1
    if dec > 0:
        raise ValueError('Requested size is too small for this value')
    return bin_str

###############################################################################
def dec2natbin(dec, n):
    r"""
    This function computes the natural binary value of a positive integer.
    
    Parameters
    ----------
    dec : number
        The decimal input integer to be converted
    n : number
        The number of bits of the output string

    Output
    ------
    bin_str : string
        The binary value in a string format

    Raises
    ------
    ValueError
        When n is too small to do the conversion

    See Also
    --------
    num : string to number conversion

    """
    bin_str = ''
    if dec >= 0:
        while n > 0:
            bin_str = str(dec % 2) + bin_str
            dec >>= 1
            n += -1
    if dec > 0:
        raise ValueError('Requested size is too small for this value')
    return bin_str

###############################################################################

