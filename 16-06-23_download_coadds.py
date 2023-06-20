#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:25:12 2023

16-06-23_download_coadds.py

GOAL:
    Use the coadd urls to download coadd images
    
    
INSTRUCTIONS:
    go to terminal
    ssh into one of the nodes (ie Summer02)
    conda activate astroconda
    python <filename>
    
    then, it should run

@author: dzakaria
"""

#%% imports 

import os
from astropy.table import QTable
import time

#%% cd into the correct directory

# os.system('cd')
# os.system('cd dzfiles/')


#%% define directory 

directory ='/users/dzakaria/dzfiles/'


#%% read in table with urls

epochs_tab=QTable.read('{path}16-06-23_epochs_table_url_size0_1.csv'.format(path=directory))

# epochs_tab=QTable.read('16-06-23_epochs_table_url_size0.1.csv')

#%%

os.system("echo coadd download beginning")



#%% loop through the table and use urls to download coadds

# for row in range(len(epochs_tab)): 
for row in range(125): 
       
        
    name = epochs_tab[row]['obj_name']
    # remove slashes
    name = name.replace('/','_')
    name = name.replace(' ', '')
    epoch = float(epochs_tab[row]['obj_epoch'])

    
    
    
    file_b1 = '{name}_ep{epoch}_b{band}.xml'.format(name=name, epoch=epoch, band=1)  
    file_b2 = '{name}_ep{epoch}_b{band}.xml'.format(name=name, epoch=epoch, band=2)  
    
#    print(file_b1.replace('_b1', ''))
        
    # if the file exists, don't redo it
    
    if os.path.exists(file_b1) and os.path.getsize(file_b1) > 1024:
        pass
    else:
    
        # every 40 processes (20 rows), have the loop sleep for 15 minutes
        # so we don't overload the neowise server with requests
        sleep_min = 15
        sleep_sec = 60*sleep_min
        # if (row) % 20 == 0:
        #     time.sleep(sleep_sec)
        
        url_b1 = epochs_tab[row]['url_band1']
        os.system('wget -O {file} \"{url}\" &'.format(file=file_b1, url=url_b1))
    
    if os.path.exists(file_b2) and os.path.getsize(file_b2) > 1024:
        pass
    else:
        url_b2 = epochs_tab[row]['url_band2']
        os.system('wget -O {file} \"{url}\" &'.format(file=file_b2, url=url_b2))
    
    
#%%

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    