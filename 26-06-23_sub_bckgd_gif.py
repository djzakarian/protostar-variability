#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:10:25 2023

@author: dzakaria
"""


#%% imports
import os
import numpy as np
from astropy.io import fits
from photutils.detection import DAOStarFinder
from photutils.aperture import aperture_photometry, CircularAperture
from astropy.visualization import ImageNormalize
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image
from astropy.wcs import WCS
import matplotlib
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils import CircularAperture
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from matplotlib.patches import Circle

# matplotlib.use('TkAgg')



#%% directory
obj_name = 'B335'
directory = f'/users/dzakaria/DATA/dzfiles/coadds-0_0667/{obj_name}/'
fits_dir_b1 = directory + '/mosaic-int_fits/b1/'
file_list_b1 = [file for file in os.listdir(fits_dir_b1) if file.endswith('.fits')]
fits_dir_b2 = directory + '/mosaic-int_fits/b2/'
file_list_b2 = [file for file in os.listdir(fits_dir_b2) if file.endswith('.fits')]

#%% 

frames = []

fig, ax = plt.subplots()

file_lists = [file_list_b1, file_list_b2]

# for file_list in file_lists:
#     for file in file_list:
#         path = directory + file

# Open the FITS file and access the data and header

# initialize first_file to track if the background has been determined
# first_file = True
for file in file_list_b1:
    
        
      
    path= directory + 'mosaic-int_fits/b1/' + file
    fig, ax = plt.subplots()
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header   
        
    center_ra = header['CRVAL1']
    center_dec = header['CRVAL2']
    pixel_scale_x = header['CDELT1']*u.deg # in deg/pixel
    pixel_scale_y = header['CDELT2'] *u.deg
    
    wcs = WCS(hdul[0].header)
    
    # Convert the pixel scale to arcseconds
    # pixel_scale_x = pixel_scale_x.to(u.arcmin).value
    pixel_scale = pixel_scale_y.to(u.arcmin).value
    
    
    
    
    if first_file == True:
        
        # Define the radius of the aperture in pixels
        ap_rad_arcmin = 0.5
        aperture_radius = ap_rad_arcmin / pixel_scale
        
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
        
    best_pixel_position = wcs.world_to_pixel_values(best_wcs_position.ra, best_wcs_position.dec)
    best_center_x = best_pixel_position[0]
    best_center_y = best_pixel_position[1]
        
    # Create the final aperture at the optimal center
    aperture = CircularAperture((best_center_x, best_center_y), aperture_radius)
    
    # change first_file to false for the next image
    first_file = False
    
    # Perform aperture photometry for the optimal aperture
    phot_table = aperture_photometry(data, aperture)
    total_ap_counts = phot_table['aperture_sum'][0]
    
    # Subtract the median background from the entire image
    data -= total_ap_counts
    
    
    # scaled_data = np.sqrt(data)
    scaled_data=data
    lower_percentile  = 0
    upper_percentile = 98
    lower_value=np.percentile(scaled_data, lower_percentile)
    upper_value=np.percentile(scaled_data, upper_percentile)
    norm = Normalize(lower_value, upper_value)
    
    # ax.cla()
    # ax.axis('off')
    
    
    
    fig.clear()
    ax = fig.add_subplot(111)
    ax.axis('off')
    
    # Display the selected region
    # fig, ax = plt.subplots()
    
    # color maps I like:
            # cubehelix, turbo, magma, inferno, gray
    im = ax.imshow(scaled_data, cmap='cubehelix',origin='lower', norm=norm)
    # aperture.plot(color='red', lw=1.5, alpha=0.7)
    aperture.plot(color='red', lw=1.5, alpha=0.7)
    
    
    
    ax.set_xlim(0, data.shape[1])
    ax.set_ylim(0, data.shape[0])
    
    canvas = FigureCanvasTkAgg(fig, master=tk.Tk())
    canvas.draw()
    canvas.get_tk_widget().pack()
    
    image = Image.frombytes('RGB', canvas.get_width_height(), fig.canvas.tostring_rgb())
    frames.append(image)
    
    
    plt.show()


#%%    

    
    
frames[0].save('/users/dzakaria/DATA/dzfiles/animations/L1527_IRAS04368+2557_b1_cubehelix.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)
    

           
#%%







