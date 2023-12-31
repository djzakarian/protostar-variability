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


#%% align images -- images is really a list of paths for fits files to be aligned


def align_images(images):
    # Read the first image to establish the reference frame
    ref_hdul = fits.open(images[0])
    reference_data = ref_hdul[0].data
    reference_header = ref_hdul[0].header

    # Initialize a list to store the aligned images
    aligned_images = [CCDData(reference_data, unit='adu', header=reference_header)]

    # Loop through the remaining images and align them to the reference frame
    for image_file in images[1:]:
        hdul = fits.open(image_file)
        data = hdul[0].data
        header = hdul[0].header

        # Perform image alignment using imreg_dft
        result = imreg_dft.imreg.similarity(reference_data, data)

        # Retrieve the aligned image and transformation parameters
        aligned_data = result['timg']

        # Create a new CCDData object with the aligned data and updated header
        aligned_image = CCDData(aligned_data, unit='adu', header=header)

        # Add the aligned image to the list
        aligned_images.append(aligned_image)

        # Close FITS file
        hdul.close()

    ref_hdul.close()

    return aligned_images




#%% make_image_user_input

def make_image_user_input(data, aperture, cmap = 'gray', plot_aperture = False,
               first_file = True, lower_percentile=0, upper_percentile=99):
    
    fig, ax = plt.subplots()
    
    
    # scaled_data = np.sqrt(data)
    scaled_data=data
    
    fig.clear()
    ax = fig.add_subplot(111)
    ax.axis('off')
    
    

    
    if first_file == True:
        
        while True:
            lower_percentile = float(input("Enter the lower percentile (0-100): "))
            upper_percentile = float(input("Enter the upper percentile (0-100): "))
            lower_value = np.percentile(scaled_data, lower_percentile)
            upper_value = np.percentile(scaled_data, upper_percentile)
            norm = Normalize(lower_value, upper_value)
    
            fig.clear()
            ax = fig.add_subplot(111)
            ax.axis('off')
    
            im = ax.imshow(scaled_data, cmap=cmap, origin='lower', norm=norm)
    
            if plot_aperture==True:
                aperture.plot(color='red', lw=1.5, alpha=0.7)
    
            ax.set_xlim(0, data.shape[1])
            ax.set_ylim(0, data.shape[0])
    
            canvas = FigureCanvasTkAgg(fig, master=tk.Tk())
            canvas.draw()
            canvas.get_tk_widget().pack()
    
            image = Image.frombytes('RGB', canvas.get_width_height(), fig.canvas.tostring_rgb())
    
            plt.close(fig)
    
            plt.imshow(image)
            plt.show()
    
            choice = input("Are you happy with the scaling? (yes/no): ")
            if choice.lower() == "yes":
                break
            
    else:
        
       
        lower_value = np.percentile(scaled_data, lower_percentile)
        upper_value = np.percentile(scaled_data, upper_percentile)
        norm = Normalize(lower_value, upper_value)

        fig.clear()
        ax = fig.add_subplot(111)
        ax.axis('off')

        im = ax.imshow(scaled_data, cmap=cmap, origin='lower', norm=norm)

        if plot_aperture==True:
            aperture.plot(color='red', lw=1.5, alpha=0.7)

        ax.set_xlim(0, data.shape[1])
        ax.set_ylim(0, data.shape[0])

        canvas = FigureCanvasTkAgg(fig, master=tk.Tk())
        canvas.draw()
        canvas.get_tk_widget().pack()

        image = Image.frombytes('RGB', canvas.get_width_height(), fig.canvas.tostring_rgb())

        plt.close(fig)

        plt.imshow(image)
        plt.show()
    
    return image, lower_percentile, upper_percentile



#%% get_background_aperture

def get_background_aperture(data, wcs, pixel_scale, radius):
    # radius units: arcmin
    
    aperture_radius = radius / pixel_scale
    
    # Calculate the minimum distance from the edge of the image in pixels
    min_distance = aperture_radius
    
    # Define the boundaries for searching centers
    min_boundary = int(min_distance )
    max_boundary_x = data.shape[0] - min_boundary
    max_boundary_y = data.shape[1] - min_boundary
    
    # Initialize variables to store the minimum median counts and corresponding center coordinates
    min_median_counts = np.inf
    best_center_x, best_center_y = 0, 0
    
    # Iterate over possible center coordinates within the boundaries
    for center_x in range(min_boundary, max_boundary_x):
        for center_y in range(min_boundary, max_boundary_y):
            # Create an aperture centered at (center_x, center_y)
            aperture = CircularAperture((center_x, center_y), aperture_radius)
    
            # Perform aperture photometry to calculate the total counts within the aperture
            phot_table = aperture_photometry(data, aperture)
            total_counts = phot_table['aperture_sum'][0]
    
            # Check if the current total counts are lower than the previous minimum
            if total_counts < min_median_counts:
                min_median_counts = total_counts
                best_center_x, best_center_y = center_x, center_y
    
    best_wcs_position = wcs.pixel_to_world(best_center_x, best_center_y) # convert to real world coords
    
    return best_wcs_position, aperture_radius
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





