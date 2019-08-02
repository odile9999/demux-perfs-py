# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
"""
    General_tools module
    ====================

    Developped by: L. Ravera

    Project: Athena X-IFU / DRE-DEMUX

    General purpose tools needed for the DRE data processing

    Routine listing
    ===============
    get_csv()
    print_conf()
    get_session_info()
    print_session_info()

    """

# -----------------------------------------------------------------------
# Imports
import os, csv
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# -----------------------------------------------------------------------
def checkdir(dirname):
    r"""
        This function checks if a directory exists. If not it creates it.

        Parameters:
        -----------
        dirname: String
        Name of the directory to be verified / created.

        Returns
        -------
        Nothing

        """
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    return()

# -----------------------------------------------------------------------
def get_csv(filename):
    r"""
        This function reads a dictionnary from a csv file.

        Parameters:
        -----------
        filename: string
        The name of the csv file

        Returns
        -------
        dictionnary: dictionnary

        """

    dictionnary={}

    if not os.path.exists(filename):
        print("File "+filename+" not found.")
    else:
        with open(filename, newline='') as csvfile:
            dict_reader = csv.reader(csvfile, delimiter=';', quotechar='|')
            for row in dict_reader:
                try:    # for numbers
                    dictionnary[row[0]]=float(row[1].replace(',','.'))
                except: # for strings
                    dictionnary[row[0]]=row[1].replace(',','.')
    return(dictionnary)

# -----------------------------------------------------------------------
def print_conf(config):
    r"""
        This function prints the parameters of a configuration set.
        (path, ...).

        Parameters:
        ----------
        config: Dictionnary

        Returns
        -------
        Nothing

        """

    print('The configuration parameters are the following:')
    for key in config.keys():
        print(key, ': ', config[key])    
    return()

# -----------------------------------------------------------------------

def get_session_info(fulldirname):
    r"""
        This function reads the session informations from a text file.

        Parameters:
        -----------
        fulldirname: string
        Name of the directory containing the text file

        Returns
        -------
        dictionnary: dictionnary

        """

    dict_info={}
    filename=os.path.join(fulldirname, 'session_informations.txt')

    if not os.path.exists(filename):
        print("Session information file not found.")
    else:
        fich = open(filename, "r")
        for line in fich:
            if '=' in line: # this is a parameter definition line
                dict_info[line.split('=')[0]]=line.split('=')[1]

    return(dict_info)

# -----------------------------------------------------------------------
def print_session_info(session_info):
    r"""
        This function prints the content of the session information dictionnary.

        Parameters:
        ----------
        session_info: Dictionnary

        Returns
        -------
        Nothing

        """

    print('The session informations are the following:')
    for key in session_info.keys():
        print(key, ': ', session_info[key])
    return()
    
# -----------------------------------------------------------------------

def non_empty_lines(table):
    r"""
        This function looks for empty lines in an 2 dimensionnal array.
    """
    n_lines = len(table[0,:])
    non_empty_lines = np.ones((n_lines), dtype=bool)
    for line in range(n_lines):
        if np.abs(table[:,line]).max()==0:
            non_empty_lines[line]=False
    return(non_empty_lines)

# ---------------------------------------------------------------------------

