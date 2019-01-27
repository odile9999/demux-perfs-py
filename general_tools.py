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

def get_hk(fulldirname, config):
    r"""
        This function gets the hk of the demux prototype from a csv file.

        Parameters:
        -----------
        fulldirname : string
        The name of the dump file (with the path)

        config : dictionnary
        Contains path and constants definitions

        Returns
        -------
        hk : dictionnary

        """

    hk={}
    hkdirname = os.path.join(fulldirname, os.path.normcase(config['dir_hk']))
    hkfilename = [f for f in os.listdir(hkdirname) \
            if os.path.isfile(os.path.join(hkdirname, f)) \
            and f[-4:]=='.csv']
    if len(hkfilename) == 0:
        print("Hk file not found.")
        hk = 0
    else:
        hkfullfilename=os.path.join(hkdirname, hkfilename[0])
        with open(hkfullfilename, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=';', quotechar='|')

            for row in reader:
                n = len(row)
                if reader.line_num == 1:
                    keys = row
                else:
                    if reader.line_num == 2:
                        for i in range(n):
                            hk[keys[i]]=np.array([row[i]])
                    else:
                        for i in range(n):
                            hk[keys[i]] = np.append(hk[keys[i]], row[i])
    return(hk)

# -----------------------------------------------------------------------

def print_hk(hk):
    r"""
        This function prints the hk of the DRE-DEMUX prototype.

        Parameters
        ----------
        hk: Dictionnary

        Returns
        -------
        Nothing

        """

    print('The hk values are the following:')
    for key in hk.keys():
        print(key, ': ', hk[key])
    
    return()
# -----------------------------------------------------------------------

def plot_hk(hk, fulldirname, config):
    r"""
        This function plots the hk of the DRE-DEMUX prototype.

        Parameters
        ----------
        hk: Dictionnary

        Returns
        -------
        Nothing

        """

    plotdirname = os.path.join(fulldirname, config['dir_plots'])
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)
    pltfilename = os.path.join(plotdirname, 'PLOT_HK.png')

    # detection of valid hks
    n_valid_hk = 0
    for key in hk.keys():
        if key[0:8]!='DRE_Hks_' and key[1:9]!='DRE_Hks_':
            n_valid_hk = n_valid_hk+1 

    n_valid_hk = n_valid_hk-1 # Date does not count as an HK

    n_cols = 3
    n_lines = n_valid_hk // n_cols
    if n_valid_hk % n_cols != 0:
        n_lines = n_lines+1

    deltatime = np.array([(datetxt_to_date(hk['Date'][0])-datetxt_to_date(hk['Date'][0])).total_seconds()])
    for i in range(len(hk['Date'])-1):
        deltatime = np.append(deltatime, (datetxt_to_date(hk['Date'][i+1])-datetxt_to_date(hk['Date'][0])).total_seconds())

    fig = plt.figure(figsize=(10, 18))
    ihk=1
    for key in hk.keys():
        if key !='Date' and ihk <= n_valid_hk:
            ax = fig.add_subplot(n_lines, n_cols, ihk)
            ax.plot(deltatime, hk[key])
            ax.grid(color='k', linestyle=':', linewidth=0.5)
            ax.set_title(key)
            ax.set_xlabel('time (s)')
            ihk=ihk+1
    fig.tight_layout()
    plt.savefig(pltfilename, bbox_inches='tight')
   
    return()
# -----------------------------------------------------------------------

def datetxt_to_date(datetxt):
    r"""
        This function converts the date string as used by the GSE into a python date.
        """
    date_test = datetime(int(datetxt[:4]),   \
                        int(datetxt[4:6]),   \
                        int(datetxt[6:8]),   \
                        int(datetxt[9:11]),  \
                        int(datetxt[11:13]), \
                        int(datetxt[13:15])) 
    return(date_test)

# -----------------------------------------------------------------------
