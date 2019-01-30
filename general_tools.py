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
    get_conf_csv()
    print_conf()

    """

# -----------------------------------------------------------------------
# Imports
import os, csv
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# -----------------------------------------------------------------------
def get_conf():
    r"""
        This function reads a dictionnary from a csv file.

        Parameters:

        filename: string
        The name of the csv file

        Returns
        -------
        dictionnary: dictionnary

        """

    dictionnary={}
    filename='demux_tools_cfg.csv'

    if not os.path.exists(filename):
        print("Configuration file not found.")
        print("It is required to define the path.")
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

        Parameters
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
