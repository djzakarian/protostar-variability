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

from astropy import wcs
from photutils import SkyCircularAperture


#%% make a gif for band 1


directory = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/L1527_IRAS04368+2557/mosaic-int_fits/b1/'
file_list = [file for file in os.listdir(directory) if file.endswith('.fits')]
           

frames = []

fig, ax = plt.subplots()

for file in file_list:
    path = directory + file
    hdulist = fits.open(path)
    data = hdulist[0].data
          
    # Get the image dimensions in pixels
    image_width = data.shape[1]
    image_height = data.shape[0]
    
    # Get the WCS information from the FITS header
    image_wcs = wcs.WCS(hdulist[0].header)
    
    # Get the central coordinates from the WCS
    ra_center, dec_center = image_wcs.all_pix2world(image_width / 2, image_height / 2, 0)
    
    # Create a SkyCoord object for the center region
    center_coords = SkyCoord(ra=ra_center, dec=dec_center, unit='deg')
    
    # Define the size of the circular region in arcseconds
    region_radius_arcsec = 4 * u.arcsec
    
    # Get the pixel scale from the FITS header
    pixel_scale_x = hdulist[0].header['CDELT1'] 
    pixel_scale_y = hdulist[0].header['CDELT2'] 
    
    # Calculate the conversion factor from arcseconds to pixels
    pixel_scale = np.mean([pixel_scale_x, pixel_scale_y]) 
    
    # Convert the region radius from arcsec to pixels
    region_radius_pixel = (region_radius_arcsec / pixel_scale)
    
    # Create a circular aperture object for the region
    region_aperture = SkyCircularAperture(center_coords, r=region_radius_pixel)

    # Perform aperture photometry to get the median counts within the circular aperture
    phot_table = aperture_photometry(data, region_aperture, method='exact')
    median_counts = phot_table['aperture_sum'][0] / region_aperture.area()



    
    scaled_data = np.sqrt(data)
    # scale_interval = PercentileInterval(0,99)
    lower_percentile  = 0
    upper_percentile = 99.25
    lower_value=np.percentile(scaled_data, lower_percentile)
    upper_value=np.percentile(scaled_data, upper_percentile)
    norm = Normalize(lower_value, upper_value)
    
    ax.cla()
    ax.axis('off')
    
    
    im=ax.imshow(scaled_data, cmap='gray',origin='lower', norm=norm)
    # Plot the circular aperture
    region_aperture.plot(color='red')
    
    fig.canvas.draw()
    image = Image.frombytes('RGB', fig.canvas.get_width_height(), fig.canvas.tostring_rgb())
    frames.append(image)
    
    # plt.imshow(scaled_data, cmap='gray', origin='lower', norm=norm)
    # plt.colorbar()
    hdulist.close()
    plt.show()
    
    
frames[0].save('animation_b1.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)


#%% make a gif for band 2

directory = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/L1527_IRAS04368+2557/mosaic-int_fits/b2/'
file_list = [file for file in os.listdir(directory) if file.endswith('.fits')]
           

frames = []

fig, ax = plt.subplots()

for file in file_list:
    path = directory + file
    hdulist = fits.open(path)
    data = hdulist[0].data
   


    
    
    scaled_data = np.sqrt(data)
    # scale_interval = PercentileInterval(0,99)
    lower_percentile  = 0
    upper_percentile = 99.25
    lower_value=np.percentile(scaled_data, lower_percentile)
    upper_value=np.percentile(scaled_data, upper_percentile)
    norm = Normalize(lower_value, upper_value)
    
    ax.cla()
    ax.axis('off')
    
    
    im=ax.imshow(scaled_data, cmap='gray',origin='lower', norm=norm)
    
    fig.canvas.draw()
    image = Image.frombytes('RGB', fig.canvas.get_width_height(), fig.canvas.tostring_rgb())
    frames.append(image)
    
    # plt.imshow(scaled_data, cmap='gray', origin='lower', norm=norm)
    # plt.colorbar()
    hdulist.close()
    plt.show()
    
    
frames[0].save('animation_b2.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)