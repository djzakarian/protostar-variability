#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 09:48:44 2023
NOTE: for some reason, this code moves all renamed files to the main directory....
I don't know why it doesn't save back to the same path


@author: dzakaria
"""

# for all files in a given directory:
    
# rename to get rid of the .xml in the middle

# put them in the correct b1 or b2 folder

""" 
# if you accidentally remove .xml when you shouldn't
import os
directory = os.getcwd()
suffix = '.xml'
for filename in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, filename)):
            new_filename = filename + suffix
            original_path = os.path.join(directory, filename)
            new_path = os.path.join(directory, new_filename)
            os.rename(original_path, new_path)

"""


# import os
# import re

# coadd_directory ='/users/dzakaria/dzfiles/coadds-0_0667'

# paths=[]
# # paths.append('/users/dzakaria/dzfiles')
# search_str = '.xml'
# pattern = re.compile(r"\.xml")

# # list of all directories:
    
# obj_names = ['B335','BHR7_IRAS08124-3422','BHR71','CB17','CB230','CB244','CB6',
#              'CB68','CG30','Ced110IRS4','CepheusE','DC303.8-14.2','HH111MMS',
#              'HH270VLA1','HH46_47','IRAS03282+3035','IRAS03292+3039','IRAS04166+2706',
#              'IRAS04169+2702','IRAS04302+2247','IRAS04325+2402','IRAS05295+1247',
#              'IRAS05329-0505','IRAS09449-5052','IRAS11072-7727','IRAS15398-3359',
#              'IRAS16253-2429','IRAS02086','L1152','L1157','L1165','L1251A','L1251B',
#              'L1251C','L1448IRS2','L1448IRS3','L1448-mm','L1489IRS','L1521F',
#              'L1527_IRAS04368+2557','L1551IRS5','L1551NE','L1616MMS1A','L1634',
#              'L483','L723_IRAS19156+1906','L778','RCrAIRAS32','Serpens1','SerpensMMS3']
    
# new_directory_list = ['results_html', 'framesused_tbl', 'mosaic-int_fits',
#                       'mosaic-unc_fits', 'mosaic-cov_fits', 'mosaic-std_fits',
#                       'mosaic-int_jpg', 'mosaic-unc_jpg', 'mosaic-cov_jpg', 
#                       'mosaic-std_jpg']
# bands = ['b1', 'b2']

# for obj in obj_names:
#     for sub_dir in new_directory_list:
#         path = ('{coadd_directory}/{obj}/{sub_dir}'.format(coadd_directory = coadd_directory, obj=obj, sub_dir=sub_dir))
#         paths.append(path)
#         for band in bands:
#             path = ('{coadd_directory}/{obj}/{sub_dir}/{band}'.format(coadd_directory=coadd_directory, obj=obj, sub_dir=sub_dir, band=band))
#             paths.append(path)


    
# for path in paths:    
#     for filename in os.listdir(path):
#         if search_str in filename:
#             new_filename = re.sub(pattern, "", filename)
#             original_path = os.path.join(path, filename)
#             new_path = os.path.join(path, new_filename)
#             os.rename(original_path, new_path)
        


