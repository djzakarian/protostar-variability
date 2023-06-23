#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 12:56:44 2023

22-06-23_download-fits-from-xml.py

GOAL: 
    Read the .xml files with download urls that allow us to download .fits files 
    of the coadds.
    
INSTRUCTIONS: RUN THIS FILE IN COMMAND LINE TO DOWNLOAD THE FILES USING WGET
                1. cd into the directory where all of the xml files are housed 
                2. conda activate astroconda
                3. python {this python code path}

@author: dzakaria
"""
#%% Make sure all subdirectories exist

# INSTRUCTIONS: Open command prompt, 
#               cd into the directory in which the subdirectories should be placed
#               type in 'python' and click Enter.
#               Copy and paste the following into the terminal

"""
# THE FOLLOWING CODE NEEDS TO BE RUN IN COMMAND LINE... 
# first: make new folders in the directory to sort the type of download 
# (results_html, framesused_tbl, mosaic-int_fits, mosaic-unc_fits, mosaic-cov_fits,
#  mosaic-std_fits, mosaic-int_jpg, mosaic-unc_jpg, mosaic-cov_jpg, mosaic-std_jpg)



new_directory_list = ['results_html', 'framesused_tbl', 'mosaic-int_fits',
                      'mosaic-unc_fits', 'mosaic-cov_fits', 'mosaic-std_fits',
                      'mosaic-int_jpg', 'mosaic-unc_jpg', 'mosaic-cov_jpg', 
                      'mosaic-std_jpg']
for new_dir in new_directory_list:
    if not os.path.exists(new_dir):
        # create the directory
        os.makedirs(new_dir)
        print('{new_dir} created sucessfully'.format(new_dir=new_dir))
    else: 
        print('{new_dir} already exists'.format(new_dir=new_dir))
"""

#%% imports

import xml.etree.ElementTree as ET
import os
from astropy.table import QTable, Column
from datetime import datetime


#%% define directory and date

directory ='/users/dzakaria/dzfiles/coadds-0_0667/'

current_datetime =  datetime.now()
date = current_datetime.strftime('%Y-%m-%d_%H:%M:%S') # update this so the log is correct


#%% SAMPLE: Read in one file and view the .xml file contents (can be left commented out)

# # make a function that prints all of the info in the .xml file 

# name = 'L1527_IRAS04368+2557_size0_0667_ep0_b1.xml'

# tree = ET.parse("{directory}{name}".format(directory=directory, name=name))
# root = tree.getroot()


# def print_elements(element, indent=''): 
#     print(f'{indent}Tag: {element.tag}')
#     print(f'{indent}Attributes: {element.attrib}')
#     if element.text:
#         print(f'{indent}Text: {element.text.strip()}')
#         for child in element:
#             print_elements(child, indent + ' ')
            
# print_elements(root)           




# # Specify the path to your XML file
# xml_file = f'{directory}{name}'

# # Specify the directory where you want to save the downloaded files
# save_directory = 'directory'

# # Parse the XML file
# tree = ET.parse(xml_file)
# root = tree.getroot()

# # Iterate over the XML elements containing the URLs
# for url_element in root.iter():
#     try:
            
#         if url_element.tag not in relevant_tags:
#             continue
        
#         url = url_element.text
#         file_name = url.split('/')[-1]  # Extract the file name from the URL
        
#         # checkpoint
#         print(f'url: {url}')
        
#         # Build the full path to save the file
#         file_path = os.path.join(save_directory, file_name)
    
#         # Download the file
#         os.system(f'wget -O "{file_path}" "{url}"')
    
#     except:
#         print('error')


#%% check which xml files have already been opened 


# make a list of all .xml files in the directory
files = os.listdir(directory)
xml_files = [file for file in files if file.endswith('.xml')]

# open and read the log of all processed xml files

# cross compare the lists, and remove any processed xml files from the xml_files list 


#%% 

log = QTable()

# how to choose acceptable dtype:
# **********************************
file_name_col=Column(xml_files, name='xml_filename')
# file_name_col.dtype='<U67'

log.add_column(file_name_col)
log['xml_filename']=xml_files
log['err_counter']=0

log.remove_rows(slice(None))


#%% 

relevant_tags=['resultHtml', 'framesused', 'images', 'fits', 'jpg', 'fits', 'jpg', 'fits', 'jpg', 'fits', 'jpg']


#%% loop through unprocessed xml files and download relevant files
rownum_counter=0
for xml_filename in xml_files:
    
    directory=directory # just a reminder that it's defined above
    
    # initialize err counter for each file
    err = 0
        
    # print rownum counter
    print(f'rownum: {rownum_counter}')
    
    
    
    # Specify the path to your XML file
    xml_file = '{directory}{filename}'.format(directory=directory, filename=xml_filename)
    
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
        file_name = xml_filename + url.split('/')[-1]  
        
        # remove the \n character at the end of file_name
        file_name = file_name.replace('\n','').strip()
        
        
        # checkpoint
        print(f'file name: {file_name}')
        print(f'url: {url}')
    
        # define save_directory depending on the file being downloaded
        if file_name.endswith('result.html'):
            save_path = directory + 'results_html/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('query_used.tbl'):
            save_path = directory + 'framesused_tbl/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-int.fits'):
            save_path = directory + 'mosaic-int_fits/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-int.jpg'):
            save_path = directory + 'mosaic-int_jpg/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-unc.fits'):
            save_path = directory + 'mosaic-unc_fits/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-unc.jpg'):
            save_path = directory + 'mosaic-unc_jpg/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-cov.fits'):
            save_path = directory + 'mosaic-cov_fits/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-cov.jpg'):
            save_path = directory + 'mosaic-cov_jpg/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-std.fits'):
            save_path = directory + 'mosaic-std_fits/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        elif file_name.endswith('mosaic-std.jpg'):
            save_path = directory + 'mosaic-std_jpg/'
            
            # Build the full path to save the file
            file_path = os.path.join(save_path, file_name)
        
            # Download the file
            os.system(f'wget -O "{file_path}" "{url}"')
            
        else: # if none of those worked, record the object in the log
            err +=1
    
    # update log
    log.add_row()
    log['xml_filename'][rownum_counter] = xml_filename
    log['err_counter'][rownum_counter] = err
    
    # save astropy table
    log.write(f'{directory}log_{date}.csv',  format='csv', overwrite=True)
    
    # update rownum counter
    rownum_counter+=1
    
    break #debugging



        

#%%

