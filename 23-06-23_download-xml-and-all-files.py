#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 16:07:07 2023

GOALS: Putting everything together:
    make a single code that goes through the coadd urls,
    downloads the relevant xml files,
    open the xml file, and download all files 
    
INSTRUCTIONS: 
    This file needs to be run on command line:
        1. cd into the directory where all of the xml files are housed (coadd directory) (IMPORTANT)
        2. conda activate astroconda
        3. python {this python code path}

@author: dzakaria
"""

#%% imports

import xml.etree.ElementTree as ET
import os
from astropy.table import QTable, Column
from datetime import datetime
import time
import numpy as np
import psutil


#%% define directory and establish date (for naming the log)

directory ='/users/dzakaria/DATA/dzfiles/'
coadd_directory ='/users/dzakaria/DATA/dzfiles/coadds-0_0667/'

current_datetime =  datetime.now()
date = current_datetime.strftime('%Y-%m-%d_%H:%M:%S') # update this so the log is correct

#%% read in the table with all of the coadd urls are held

size ='0_0667'

epochs_tab=QTable.read('{path}29-06-23_epochs_table_url_size{size}.csv'.format(path=directory, size=size))


#%% make directories for each object

# make list of target names
obj_names = np.unique(epochs_tab['obj_name'])

# remove slashes and spaces
for i in range(len(obj_names)):
    obj_names[i] = obj_names[i].replace('/','_')
    obj_names[i] = obj_names[i].replace(' ', '')



for new_dir in obj_names:
    if not os.path.exists(new_dir):
        # create the directory
        os.makedirs(new_dir)
        print(f'{new_dir} created sucessfully')
    else: 
        print(f'{new_dir} already exists')
        
        
        
#%% make new subdirectories for each file type within each object directory

new_directory_list = ['results_html', 'framesused_tbl', 'mosaic-int_fits',
                      'mosaic-unc_fits', 'mosaic-cov_fits', 'mosaic-std_fits',
                      'mosaic-int_jpg', 'mosaic-unc_jpg', 'mosaic-cov_jpg', 
                      'mosaic-std_jpg']

for obj in obj_names:
    for sub_dir in new_directory_list:
        if not os.path.exists(f'{obj}/{sub_dir}'):
            os.makedirs(f'{obj}/{sub_dir}')
            # print(f'{sub_dir} created sucessfully')
        # else: 
            # print(f'{sub_dir} already exists')
 
#%% now put sub-subdirectories that sort each file by band 
bands = ['b1','b2']
for obj in obj_names:
    for sub_dir in new_directory_list:
        for band in bands:
            if not os.path.exists(f'{obj}/{sub_dir}/{band}'):
                os.makedirs(f'{obj}/{sub_dir}/{band}')


#%% 

relevant_tags=['resultHtml', 'framesused', 'images', 'fits', 'jpg', 'fits', 'jpg', 'fits', 'jpg', 'fits', 'jpg']



#%% loop through the epochs_tab table and use urls to download coadds and relevant files

for row in range(0, len(epochs_tab)): #note: come back to row 129 bc it wasn't working

    
    name = epochs_tab[row]['obj_name']
    # remove slashes and spaces
    name = name.replace('/','_')
    name = name.replace(' ', '')
    epoch = int(epochs_tab[row]['obj_epoch'])

    

    
    file_b1 = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch, band=1)  
    file_b2 = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch, band=2)  
    

        


    
    if os.path.exists(file_b1) and os.path.getsize(file_b1) > 1024:
        pass # check if the next file exists
    else:
        
        
        
        # checkpoint
        print(f'row: {row}')

        
        # check in 
        print(file_b1.replace('_b1', ''))

        
        url_b1 = epochs_tab[row]['url_band1']
        
        
        # start the timer
        start_time = time.time()
        os.system(f'wget -O {file_b1} -t 3 \"{url_b1}\"')
 
        if os.path.exists(file_b1) and os.path.getsize(file_b1) > 1024:
            pass # keep going in this row, download was succesful
        
        else:
             continue # move on, unsuccessful download
    
    
    if os.path.exists(file_b2) and os.path.getsize(file_b2) > 1024:
        continue # move on to next row 
    else:
        
            
        

        
        url_b2 = epochs_tab[row]['url_band2']
        
        
        # start the timer
        start_time = time.time()
        os.system(f'wget -O {file_b2} -t 3 \"{url_b2}\"')
 
        if os.path.exists(file_b2) and os.path.getsize(file_b2) > 1024:
            pass # keep going in this row, download was succesful
        
        else:
             continue # move on, unsuccessful download
    
            
        
        
    # now that the xml is downloaded, get the files downloaded from the xml
        
    xml_files = [file_b1, file_b2]
    for xml_file in xml_files:
        # checkpoint
        print(f'xml_file: {xml_file}')
        

        # Parse the XML file
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # checkpoint
        print(tree) 
        
        # Iterate over the XML elements containing the URLs
        for url_element in root.iter():
            
            # first: ignore the elements that don't have the file info
            if url_element.tag not in relevant_tags:
                continue
            
            # checkpoint
            print(f'url_element: {url_element}')
            
            url = url_element.text
            
            # remove '/n' and extra space at the end
            url = url.replace('\n', '').strip()
            
            # NAMING: first part is the og filename, 
            # then extract the file name from the URL and append to the end
            # debugging: this used to be xml_filename instead of filename... make sure this was the correct way to fix code
            file_name = xml_file + url.split('/')[-1]  
            
            # remove the \n character at the end of file_name
            file_name = file_name.replace('\n','').strip()
            
            # remove the .xml
            file_name = file_name.replace('.xml', '').strip()
            
            
            # checkpoint
            print(f'file name: {file_name}')
            print(f'url: {url}')
        
            # html doesn't save useful info
            # # define save_directory depending on the file being downloaded
            # if file_name.endswith('b1result.html'):
            #     save_path = coadd_directory + f'{name}/' + 'results_html/' + 'b1/'
                
                # # Build the full path to save the file
                # file_path = os.path.join(save_path, file_name)
            
                # # Download the file
                # os.system(f'wget -O "{file_path}" "{url}"')
                
            if file_name.endswith('W1_query_used.tbl'):
                save_path = coadd_directory + f'{name}/' + 'framesused_tbl/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-int.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-int_fits/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-int.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-int_jpg/' + 'b1/'
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-unc.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-unc_fits/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-unc.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-unc_jpg/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-cov.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-cov_fits/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-cov.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-cov_jpg/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-std.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-std_fits/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W1_mosaic-std.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-std_jpg/' + 'b1/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
    
            # html doesn't save useful info
            # elif file_name.endswith('b2result.html'):
            #     save_path = coadd_directory + f'{name}/' + 'results_html/' + 'b2/'
                
            #     # Build the full path to save the file
            #     file_path = os.path.join(save_path, file_name)
            
            #     # Download the file
            #     os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_query_used.tbl'):
                save_path = coadd_directory + f'{name}/' + 'framesused_tbl/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-int.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-int_fits/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-int.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-int_jpg/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-unc.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-unc_fits/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-unc.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-unc_jpg/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-cov.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-cov_fits/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-cov.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-cov_jpg/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-std.fits'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-std_fits/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
                
            elif file_name.endswith('W2_mosaic-std.jpg'):
                save_path = coadd_directory + f'{name}/' + 'mosaic-std_jpg/' + 'b2/'
                
                # Build the full path to save the file
                file_path = os.path.join(save_path, file_name)
            
                # Download the file
                os.system(f'wget -O "{file_path}" "{url}"')
        

