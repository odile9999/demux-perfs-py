# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------
"""
    hk_tools module
    ====================

    Developped by: L. Ravera

    Project: Athena X-IFU / DRE-DEMUX

    General purpose tools needed for the DRE data processing

    Routine listing
    ===============

    """

# -----------------------------------------------------------------------
# Imports
import os, csv
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# -----------------------------------------------------------------------

def get_hk(fulldirname, config):
    r"""
        This function gets the hk of the demux prototype from a csv file.

        Parameters:
        -----------
        fulldirname: string
        The name of the dump file (with the path)

        config: dictionnary
        Contains path and constants definitions

        Returns
        -------
        hk : dictionnary

        """

    hk=hk_lims={}
    hkdirname = os.path.join(os.path.normcase(fulldirname), os.path.normcase(config['dir_hk']))
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
                            if keys[i] == 'Date':
                                hk[keys[i]]=np.array([row[i]])
                            else:
                                hk[keys[i]]=np.array([float(row[i].replace(',','.'))])
                    else:
                        for i in range(n):
                            if keys[i] == 'Date':
                                hk[keys[i]] = np.append(hk[keys[i]], row[i])
                            else:
                                hk[keys[i]] = np.append(hk[keys[i]], float(row[i].replace(',','.')))
        hk_lims = get_hk_lims(fulldirname, config, hk)
    return(hk, hk_lims)

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

def plot_hk(hk, hk_lims, fulldirname, config):
    r"""
        This function plots the hk of the DRE-DEMUX prototype.

        Parameters
        ----------
        hk: Dictionnary
        hk values.

        hk_lims: Dictionnary
        hk limits. 

        fulldirname: string
        The name of the dump file (with the path)

        config: dictionnary
        Contains path and constants definitions

        Returns
        -------
        Nothing

        """

    plotdirname = os.path.join(os.path.normcase(fulldirname), os.path.normcase(config['dir_plots']))
    if not os.path.isdir(plotdirname):
        os.mkdir(plotdirname)
    pltfilename1 = os.path.join(plotdirname, 'PLOT_HK.png')
    pltfilename2 = os.path.join(plotdirname, 'PLOT_HK_LIMS.png')

    # detection of valid hks
    n_valid_hk = 0
    for key in hk.keys():
        if key[0:8]!='DRE_Hks_' and key[1:9]!='DRE_Hks_':
            n_valid_hk = n_valid_hk+1 

    n_valid_hk = n_valid_hk-1 # Date does not count as an HK

    n_cols = 4
    n_lines = n_valid_hk // n_cols
    if n_valid_hk % n_cols != 0:
        n_lines = n_lines+1

    deltatime = np.array([(datetxt_to_date(hk['Date'][0])-datetxt_to_date(hk['Date'][0])).total_seconds()])
    for i in range(len(hk['Date'])-1):
        deltatime = np.append(deltatime, (datetxt_to_date(hk['Date'][i+1])-datetxt_to_date(hk['Date'][0])).total_seconds())

    fig = plt.figure(figsize=(12, 18))
    ihk=1
    for key in hk.keys():
        if key !='Date' and ihk <= n_valid_hk:
            ax = fig.add_subplot(n_lines, n_cols, ihk)
            ax.plot(deltatime, hk[key], linewidth=3)
            ax.grid(color='k', linestyle=':', linewidth=0.5)
            ax.set_title(key[1:-1])
            ax.set_xlabel('time (s)')
            if hk_lims[key]['lowalert']:
                ax.plot([deltatime[0], deltatime[-1]], [hk_lims[key]['lowalertv'], hk_lims[key]['lowalertv']], '--', color='red')
            if hk_lims[key]['lowwarn']:
                ax.plot([deltatime[0], deltatime[-1]], [hk_lims[key]['lowwarnv'], hk_lims[key]['lowwarnv']], '--', color='orange')
            if hk_lims[key]['highwarn']:
                ax.plot([deltatime[0], deltatime[-1]], [hk_lims[key]['highwarnv'], hk_lims[key]['highwarnv']],'--', color='orange')
            if hk_lims[key]['highalert']:
                ax.plot([deltatime[0], deltatime[-1]], [hk_lims[key]['highalertv'], hk_lims[key]['highalertv']],'--', color='red')
            mini=min(hk[key].min(), hk_lims[key]['lowalertv'])
            maxi=max(hk[key].max(), hk_lims[key]['highalertv'])
            margin = 0.2*(maxi-mini)
            ax.set_ylim(mini-margin, maxi+margin)
            ax.grid(color='k', linestyle=':', linewidth=0.5)
            ihk=ihk+1
    fig.tight_layout()
    plt.savefig(pltfilename2, bbox_inches='tight')

    fig = plt.figure(figsize=(12, 18))
    ihk=1
    for key in hk.keys():
        if key !='Date' and ihk <= n_valid_hk:
            ax = fig.add_subplot(n_lines, n_cols, ihk)
            ax.plot(deltatime, hk[key], linewidth=3)
            ax.grid(color='k', linestyle=':', linewidth=0.5)
            ax.set_title(key[1:-1])
            ax.set_xlabel('time (s)')
            mini=hk[key].min()
            maxi=hk[key].max()
            margin = 0.2*(maxi-mini)
            ax.set_ylim(mini-margin, maxi+margin)
            ax.grid(color='k', linestyle=':', linewidth=0.5)
            ihk=ihk+1
    fig.tight_layout()
    plt.savefig(pltfilename1, bbox_inches='tight')

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

def get_hk_lims(fulldirname, config, hk):
    r"""
        This function gets the hk limits from a text file.

        Parameters
        ----------
        fulldirname : string
        The name of the dump file (with the path)

        config : dictionnary
        Contains path and constants definitions

        hk: Dictionnary
        hk values.

        Returns
        -------
        hk_lims: Dictionnary
        hk limits. 

        """

    hk_lims_dict={\
            'lowalert': False, \
            'lowalertv': 0, \
            'lowwarn': False, \
            'lowwarnv': 0, \
            'highwarn': False, \
            'highwarnv': 0, \
            'highalert': False, \
            'highalertv': 0
            }

    hk_lims_list={}
    for key in hk.keys():
        if key!='Date' and key[1:9]!="DRE_Hks_":
            hk_lims_list[key]=hk_lims_dict.copy()

    filename = 'parametersTF.dispatcher'

    hkdirname = os.path.join(os.path.normcase(fulldirname), os.path.normcase(config['dir_hk']))
    fullfilename = os.path.join(hkdirname, filename)
    if not os.path.isfile(fullfilename):
        print(filename+" file not found.")
    else:
        for key in hk.keys():
            if key!='Date' and key[1:9]!="DRE_Hks_":
                fich = open(fullfilename, "r")
                # looking for hk informations in the file
                line = ''
                hk_real_name = key[1:-5]
                while line != 'realname='+hk_real_name+'\n':
                    line=fich.readline()
                # skipping unwanted lines
                [fich.readline() for i in range(5)]
                # getting hk limits
                hk_lims_list[key]['lowalert']=fich.readline().split("=")[1]=='true\n'
                hk_lims_list[key]['lowalertv']=float(fich.readline().split('=')[1].replace(',','.'))
                hk_lims_list[key]['lowwarn']=fich.readline().split('=')[1]=='true\n'
                hk_lims_list[key]['lowwarnv']=float(fich.readline().split('=')[1].replace(',','.'))
                hk_lims_list[key]['highwarn']=fich.readline().split('=')[1]=='true\n'
                hk_lims_list[key]['highwarnv']=float(fich.readline().split('=')[1].replace(',','.'))
                hk_lims_list[key]['highalert']=fich.readline().split('=')[1]=='true\n'
                hk_lims_list[key]['highalertv']=float(fich.readline().split('=')[1].replace(',','.'))
                fich.close()

    return(hk_lims_list)
# -----------------------------------------------------------------------

def print_hk_lims(hk_lims_list):
    r"""
        This function prints the hk limits.

        Parameters
        ----------
        hk_lims: Dictionnary
        hk limits. 

        Returns
        -------
        Nothing.
        
        """

    for key1 in hk_lims_list.keys():
        if key1!='Date' and key1[1:9]!="DRE_Hks_":
            print('-------------------------------')
            print(key1)
            for key2 in hk_lims_list[key1].keys():
                print(key2, ': ', hk_lims_list[key1][key2])
# -----------------------------------------------------------------------
