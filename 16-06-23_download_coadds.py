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

for row in range(20):
    
    name = epochs_tab[row]['obj_name']
    # remove slashes
    name = name.replace('/','_')
    epoch = float(epochs_tab[row]['obj_epoch'])

    
    print('object: ', name)
    
    file_b1 = '{name}_ep{epoch}_b{band}.xml'.format(name=name, epoch=epoch, band=1)  
    file_b2 = '{name}_ep{epoch}_b{band}.xml'.format(name=name, epoch=epoch, band=2)  
    
        
    # if the file exists, don't redo it
    
    if os.path.exists(file_b1):
        continue
    else:
        url_b1 = epochs_tab[row]['url_band1']
        os.system('wget -O {file} \"{url}\" &'.format(file=file_b1, url=url_b1))
    
    if os.path.exists(file_b2):
        continue
    else:
        url_b2 = epochs_tab[row]['url_band2']
        os.system('wget -O {file} \"{url}\" &'.format(file=file_b2, url=url_b2))
    
    
#%%

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    