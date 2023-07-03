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

from photutils import SkyCircularAperture


#%% make a gif for band 1 


directory = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/L1527_IRAS04368+2557/mosaic-int_fits/b1/'
file_list = [file for file in os.listdir(directory) if file.endswith('.fits')]
           
#%% 
frames = []

fig, ax = plt.subplots()

for file in file_list:
    path = directory + file
    
    # Open the FITS file and access the data and header
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    
   
    
   

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

    plt.show()
    
    
# frames[0].save('animation_b1.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)


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
    
    
frames[0].save('/users/dzakaria/DATA/dzfiles/coadds-0_0667/L1527_IRAS04368+2557/animation_b2.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)