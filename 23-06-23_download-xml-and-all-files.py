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
from astropy.table import QTable
from astropy.coordinates import SkyCoord 
import astropy.units as u
from astropy.time import Time
from datetime import datetime
import numpy as np



#%% define directory and establish date (for naming the log)

directory ='/users/dzakaria/DATA/dzfiles/'
coadd_directory ='/users/dzakaria/DATA/dzfiles/coadds-0_0667/'


#%% read in the table with all of the coadd urls are held

size ='0_0667'

epochs_tab=QTable.read('{path}30-06-23_epochs_table_url_size{size}.csv'.format(path=directory, size=size))


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


#%% 

def sub_epoch_urls(epochs_tab_row, band, size = 0.0667):
    
    # split epoch into two and get urls
    ra = epochs_tab[row]['ra']
    dec = epochs_tab[row]['dec']
    date1_str = epochs_tab[row]['date_obs1']
    date2_str = epochs_tab[row]['date_obs2']

    
    
    # turn mjd date of middle date to a consistent format (ISO 8601 format according to chatgpt :)
    date1, date2 = epochs_tab[row]['mjd_obs1'], epochs_tab[row]['mjd_obs2']
    mid_obs = (date1 + date2) / 2
    t = Time(mid_obs, format='mjd', scale='utc')
    datemid_str = t.iso
    


    
    # formate date1 and date2 for the url
    
    # date1:
    dt1 = datetime.strptime(date1_str, '%Y-%m-%d %H:%M:%S.%f')
    date1 = dt1.strftime('%d%b%Y %H:%M:%S')
    
    dtmid = datetime.strptime(datemid_str, '%Y-%m-%d %H:%M:%S.%f') 
    datemid =  dtmid.strftime('%d%b%Y %H:%M:%S')
    
    dt2 = datetime.strptime(date2_str, '%Y-%m-%d %H:%M:%S.%f')
    date2 = dt2.strftime('%d%b%Y %H:%M:%S')
    
    coordinate = SkyCoord(ra=ra, dec=dec, unit=u.deg)
    
    coord = coordinate.to_string('hmsdms', sep=':', precision = 3)
        
    
    url_a = "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr={coord}&band={band}&sizeX={sizeX}&sizeY={sizeY}&date1={date1}&date2={date2}&mode=PI" \
            .format(coord=coord,band=band, sizeX=size, sizeY=size, date1=date1, date2=datemid)
    
    url_b = "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr={coord}&band={band}&sizeX={sizeX}&sizeY={sizeY}&date1={date1}&date2={date2}&mode=PI" \
            .format(coord=coord,band=band, sizeX=size, sizeY=size, date1=datemid, date2=date2)    
    return url_a, url_b


#%% loop through the epochs_tab table and use urls to download coadds and relevant files

for row in range(0, len(epochs_tab)): #note: come back to row 129 bc it wasn't working

    
    name = epochs_tab[row]['obj_name']
    # remove slashes and spaces
    name = name.replace('/','_')
    name = name.replace(' ', '')
    epoch = int(epochs_tab[row]['obj_epoch'])

    

    
    file_b1 = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch, band=1)  
    file_b2 = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch, band=2)  
    

        
    # initialize list of successfully downloaded xml files
    xml_files = []
    
    if os.path.exists(file_b1) and os.path.getsize(file_b1) > 1024:
        pass # check if the next file exists
        
    else: 
        
        
        # checkpoint
        print(f'row: {row}')
        # check in 
        print(file_b1.replace('_b1', ''))

        # read in url and make wget command
        url_b1 = epochs_tab[row]['url_band1']
        os.system(f'wget -O {file_b1} -t 3 \"{url_b1}\"')
 
        if os.path.exists(file_b1) and os.path.getsize(file_b1) > 1024:
            xml_files.append(file_b1)
            pass # keep going in this row, download was succesful
        
        else: 
            # epoch download was not successful
        
            
            # don't give up: try splitting the epoch in half to reduce strain on the servers
            # now there are two sub_epochs
            # the new epoch will be epoch + 0.5
            
            epoch_a = float(epoch)
            epoch_b = epoch_a + 0.5
            
            file_b1a = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch_a, band=1) 
            file_b1b = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch_b, band=1) 
            
            url_b1a, url_b1b = sub_epoch_urls(row, band=1)
          
            # attempt to download the files
            os.system(f'wget -O {file_b1a} -t 3 \"{url_b1a}\"')
            os.system(f'wget -O {file_b1b} -t 3 \"{url_b1b}\"')
            
            # now check if the new files downloaded successfully
            if (os.path.exists(file_b1a) and os.path.getsize(file_b1a) > 1024 
                and os.path.exists(file_b1b) and os.path.getsize(file_b1b) > 1024):
                xml_files.append(file_b1a)
                xml_files.append(file_b1b)
                
            else:
                continue # still didn't work, so move on 
                

    
    if os.path.exists(file_b2) and os.path.getsize(file_b2) > 1024:
        
        continue # move on to next row, these files have been downloaded
    else:
        
        # read in url and make wget command
        url_b2 = epochs_tab[row]['url_band2']
        os.system(f'wget -O {file_b2} -t 3 \"{url_b2}\"')
 
        if os.path.exists(file_b2) and os.path.getsize(file_b2) > 1024:
            xml_files.append(file_b2)
            pass # keep going in this row, download was succesful
        
        else: 
            # epoch download was not successful

            
            # don't give up: try splitting the epoch in half to reduce strain on the servers
            # now there are two sub_epochs
            # the new epoch will be epoch + 0.5
            
            epoch_a = float(epoch)
            epoch_b = epoch_a + 0.5
            
            file_b2a = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch_a, band=2) 
            file_b2b = '{name}_size{size}_ep{epoch}_b{band}.xml'.format(name=name, size=size, epoch=epoch_b, band=2) 
            
            url_b2a, url_b2b = sub_epoch_urls(row, band=2)
          
            # attempt to download the files
            os.system(f'wget -O {file_b2a} -t 3 \"{url_b2a}\"')
            os.system(f'wget -O {file_b2b} -t 3 \"{url_b2b}\"')
            
            # now check if the new files downloaded successfully
            if (os.path.exists(file_b2a) and os.path.getsize(file_b2a) > 1024 
                and os.path.exists(file_b2b) and os.path.getsize(file_b2b) > 1024):
                xml_files.append(file_b2a) 
                xml_files.append(file_b2b) 
                
            else:
                
                continue # still didn't work, so move on 
                

    
            
        
        
    # now that the xml is downloaded, get the files downloaded from the xml
        
    
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
        

