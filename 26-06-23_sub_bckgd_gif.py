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

from photutils import CircularAperture

matplotlib.use('Qt5Agg')

#%% onclick 
# global variables
x_center = None
y_center = None

def onclick(event):
    x_center = None
    y_center = None
    if event.button ==1: # left mouse button
        if x_center is None or y_center is None: # first click is center
            # get clicked coordinates
            x_center, y_center = event.xdata, event.ydata
        
        else: # second click chooses radius
            radius = np.sqrt((event.xdata-x_center)**2+(event.ydata-y_center)**2)
            
    # create the circular aperture
    aperture = CircularAperture(x_center, y_center, radius)
    
    # clear the previous plot
    plt.clf()
#%% directory
obj_name = 'L1527_IRAS04368+2557'
directory = f'/users/dzakaria/DATA/dzfiles/coadds-0_0667/{obj_name}/'
fits_dir_b1 = directory + '/mosaic-int_fits/b1/'
file_list_b1 = [file for file in os.listdir(fits_dir_b1) if file.endswith('.fits')]
fits_dir_b2 = directory + '/mosaic-int_fits/b2/'
file_list_b2 = [file for file in os.listdir(fits_dir_b2) if file.endswith('.fits')]

#%% 

frames = []

file_lists = [file_list_b1, file_list_b2]

# for file_list in file_lists:
#     for file in file_list:
#         path = directory + file

# Open the FITS file and access the data and header
for file in file_list_b1:
    path= directory + '/mosaic-int_fits/b2/' + file
    fig, ax = plt.subplots()
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header   
        
    center_ra = header['CRVAL1']
    center_dec = header['CRVAL2']
    pixel_scale = header['CDELT2']
    
    # convert pixel scale from arcsec to deg
    pixel_scale_deg = pixel_scale / 3600.0

   
    
    
    # scale the data

    scaled_data = np.sqrt(data)
    # scale_interval = PercentileInterval(0,99)
    lower_percentile  = 0
    upper_percentile = 99.25
    lower_value=np.percentile(scaled_data, lower_percentile)
    upper_value=np.percentile(scaled_data, upper_percentile)
    norm = Normalize(lower_value, upper_value)
    
    
    # plot the image
    plt.figure()
    plt.imshow(scaled_data, origin='lower', cmap='gray', norm=norm)

    
    # add labels and title
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.title('')
    
    # connect the onclick event handler
    plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    
    # display plot
    plt.show()
    
    # break

           
#%% 
frames = []

fig, ax = plt.subplots()

for file in file_list_b1:
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