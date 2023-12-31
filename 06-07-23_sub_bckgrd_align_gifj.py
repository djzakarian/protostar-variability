#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 21:03:00 2023

Generate gifs from fits files. The background aperture is automatically found
and background counts are subtracted.

@author: dzakaria
"""

#%% imports

import os
import numpy as np
from astropy.io import fits
from photutils.aperture import aperture_photometry, CircularAperture
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image
from astropy.wcs import WCS
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.wcs import WCS
from scipy.ndimage import fourier_shift
from scipy.ndimage import measurements
import numpy as np

import numpy as np
from astropy.io import fits
from astropy.nddata import CCDData
from scipy.ndimage import fourier_shift

from image_registration import chi2_shift, cross_correlation_shifts
# import cv2
# from fits_align import align
from astropy.wcs.utils import pixel_to_pixel
from reproject import reproject_exact
from scipy.ndimage import map_coordinates
import imreg_dft




#%% list of objects, bands, cmaps

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

cmaps = ['gray', 'magma', 'turbo', 'cubehelix']
cmaps=['turbo']


#%% coadd directory

coadd_directory = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/'



#%%
file_lists=[]

ap_rad_arcmin = 0.5

lower_percentile=0 # initialize normalization parameters
upper_percentile=99.5 

for name in obj_names:
    for band in bands:
        for cmap in cmaps:
            file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
            file_list = [file for file in os.listdir(file_dir) if file.endswith('.fits')]
            
            file_list_paths = []
            for file in file_list:
                file_list_paths.append(f'{file_dir}/{file}')
            
            aligned_images = align_images(file_list_paths)
            
            
            frames=[]
            first_file = True # first file in directory is used to find background aperture
            for image in aligned_images:
                
                
                fig, ax = plt.subplots()
                data = image.data
                header = image.header
                   
                center_ra = header['CRVAL1']
                center_dec = header['CRVAL2']
                pixel_scale_x = header['CDELT1']*u.deg # in deg/pixel
                pixel_scale_y = header['CDELT2'] *u.deg
                
                wcs = WCS(image[0].header)
                
                # Convert the pixel scale to arcseconds
                # pixel_scale_x = pixel_scale_x.to(u.arcmin).value
                pixel_scale = pixel_scale_y.to(u.arcmin).value
                
                
                # if loop is on the first file in a directory, 
                # then go through the image and find a suitable aperture for 
                # the background
                if first_file == True:
                    # ap_rad_arcmin initialized outside of loop
                    best_wcs_position, aperture_radius = get_background_aperture(data=data, wcs=wcs,
                                                                                pixel_scale=pixel_scale,
                                                                                radius=ap_rad_arcmin)
                    
                    
                    
                    
                # convert the pixel position to ra and dec    
                best_pixel_position = wcs.world_to_pixel_values(best_wcs_position.ra, best_wcs_position.dec)
                best_center_x = best_pixel_position[0]
                best_center_y = best_pixel_position[1]
                    
                # Create the final aperture at the optimal center
                aperture = CircularAperture((best_center_x, best_center_y), aperture_radius)
                
                
                
                # Perform aperture photometry for the optimal aperture
                phot_table = aperture_photometry(data, aperture)
                total_ap_counts = phot_table['aperture_sum'][0]
                
                # Subtract the median background from the entire image
                data -= total_ap_counts
                
                print(total_ap_counts)
                
                image, lower_percentile, upper_percentile = make_image_user_input(data=data, aperture=aperture,
                                          cmap=cmap, plot_aperture=True,
                                         first_file=first_file, lower_percentile=lower_percentile,
                                         upper_percentile=upper_percentile)
                frames.append(image)
                
                
                plt.show()
           
                # change first_file tracker to false for the next image
                first_file = False
                
            frames[0].save(f'/users/dzakaria/DATA/dzfiles/animations/{name}_{band}_{cmap}_aligned.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)   

           

#%%





