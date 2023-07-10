#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:32:30 2023

GOAL: 
    sample code using ps_var_lib


@author: dzakaria
"""

#%% list of objects, bands, cmaps

coadd_directory = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/'

obj_names = ['B335','BHR7_IRAS08124-3422','BHR71','CB17','CB230','CB244','CB6',
            'CB68','CG30','Ced110IRS4','CepheusE','DC303.8-14.2','HH111MMS',
            'HH270VLA1','HH46_47','IRAS03282+3035','IRAS03292+3039','IRAS04166+2706',
            'IRAS04169+2702','IRAS04302+2247','IRAS04325+2402','IRAS05295+1247',
            'IRAS05329-0505','IRAS09449-5052','IRAS11072-7727','IRAS15398-3359',
            'IRAS16253-2429','IRAS02086','L1152','L1157','L1165','L1251A','L1251B',
            'L1251C','L1448IRS2','L1448IRS3','L1448-mm','L1489IRS','L1521F',
            'L1527_IRAS04368+2557','L1551IRS5','L1551NE','L1616MMS1A','L1634',
            'L483','L723_IRAS19156+1906','L778','RCrAIRAS32','Serpens1','SerpensMMS3']

obj_names = ['L1527_IRAS04368+2557']

bands = ['b1', 'b2']
bands = ['b1']   
#%% code to loop through all of the objects, subtract background
    

for name in obj_names:
    for band in bands:
        file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
        # process_images(file_dir, 5, cmap = 'turbo')
        
        point_source_photometry(file_dir, lower_percentile=20, upper_percentile=99.5, fwhm=10, threshold = 7)
        
#%% 