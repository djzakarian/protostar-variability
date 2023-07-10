#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:38:39 2023

https://github.com/viveikjha/aperture_photometry/blob/master/03_aperture_photometry.ipynb

@author: dzakaria
"""

#%% imports
from astropy.io import fits
import astropy.units as u
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import ImageNormalize, LogStretch
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import numpy as np
from photutils.background import Background2D, MeanBackground
import os
from matplotlib.colors import Normalize
from astropy.wcs import WCS
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from PIL import Image
from photutils.aperture import aperture_photometry, CircularAperture
# Show plots in the notebook
# %matplotlib inline


#%

#%% list of objects, bands, cmaps

# obj_names = ['B335','BHR7_IRAS08124-3422','BHR71','CB17','CB230','CB244','CB6',
#               'CB68','CG30','Ced110IRS4','CepheusE','DC303.8-14.2','HH111MMS',
#               'HH270VLA1','HH46_47','IRAS03282+3035','IRAS03292+3039','IRAS04166+2706',
#               'IRAS04169+2702','IRAS04302+2247','IRAS04325+2402','IRAS05295+1247',
#               'IRAS05329-0505','IRAS09449-5052','IRAS11072-7727','IRAS15398-3359',
#               'IRAS16253-2429','IRAS02086','L1152','L1157','L1165','L1251A','L1251B',
#               'L1251C','L1448IRS2','L1448IRS3','L1448-mm','L1489IRS','L1521F',
#               'L1527_IRAS04368+2557','L1551IRS5','L1551NE','L1616MMS1A','L1634',
#               'L483','L723_IRAS19156+1906','L778','RCrAIRAS32','Serpens1','SerpensMMS3']
obj_names = ['L1527_IRAS04368+2557']

# bands = ['b1', 'b2']
bands = ['b1']

# cmaps = ['gray', 'magma', 'turbo', 'cubehelix']
cmaps=['magma']


#%% coadd directory

coadd_directory = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/'



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
    
#%% 

def make_image_user_input(data, aperture, cmap = 'gray', plot_aperture = True,
               first_file = True, lower_percentile=0, upper_percentile=99):
    
    fig, ax = plt.subplots()
    
    
    # scaled_data = np.sqrt(data)
    scaled_data=data
    
    
    # lower_percentile  = 0
    # upper_percentile = 100
    # lower_value=np.percentile(scaled_data, lower_percentile)
    # upper_value=np.percentile(scaled_data, upper_percentile)
    # norm = Normalize(lower_value, upper_value)

    
    fig.clear()
    ax = fig.add_subplot(111)
    ax.axis('off')
    
    
    # im = ax.imshow(scaled_data, cmap=cmap, origin='lower', norm=norm)
    
    # if plot_aperture == True:
    #     aperture.plot(color='red', lw=1.5, alpha=0.7) 

    
    
    
    # ax.set_xlim(0, data.shape[1])
    # ax.set_ylim(0, data.shape[0])
    
    # canvas = FigureCanvasTkAgg(fig, master=tk.Tk())
    # canvas.draw()
    # canvas.get_tk_widget().pack()
    
    # image = Image.frombytes('RGB', canvas.get_width_height(), fig.canvas.tostring_rgb())
    
    # plt.show()
    
    # return image
    
    
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

        # ax.set_xlim(0, data.shape[1])
        # ax.set_ylim(0, data.shape[0])

        canvas = FigureCanvasTkAgg(fig, master=tk.Tk())
        canvas.draw()
        canvas.get_tk_widget().pack()

        image = Image.frombytes('RGB', canvas.get_width_height(), fig.canvas.tostring_rgb())

        plt.close(fig)

        plt.imshow(image)
        plt.show()
    
    return image, lower_percentile, upper_percentile
#%% 

ap_rad_arcmin = 0.5 # background aperture radius in arcsc

for name in obj_names:
    for band in bands:
    
        file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
        file_list = [file for file in os.listdir(file_dir) if file.endswith('.fits')]
        
        for file in file_list:
            
            path= f'{file_dir}/{file}'
            fig, ax = plt.subplots()
            with fits.open(path) as hdul:
                data = hdul[0].data
                header = hdul[0].header 
            
            
            wcs = WCS(hdul[0].header)
            pixel_scale_y = header['CDELT2'] *u.deg
            
            # Convert the pixel scale to arcseconds
            # pixel_scale_x = pixel_scale_x.to(u.arcmin).value
            pixel_scale = pixel_scale_y.to(u.arcmin).value
            
            # define the mask
            mask = (data == 0)
            
            # store image 
            unit = u.electron/u.s
            xdf_image=CCDData(data, unit=unit, meta=header, mask=mask)
            
            
            """ DETERMINE SUITABLE BACKGROUND REGION AND SUBTRACT FROM ENTIRE IMAGE """       
            
            radius=0.5
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
        
            image, lower_percentile, upper_percentile = make_image_user_input(data=data, aperture=aperture,
                                      cmap='gray', plot_aperture=True,
                                      lower_percentile=0,
                                     upper_percentile=99.5)
            
            

            
            

            """ CREATING COMP STAR APERTURES """
            
            from photutils import find_peaks
            from photutils.centroids import centroid_2dg
        
            # calculate statistics
            mean, median, std = sigma_clipped_stats(xdf_image.data, sigma=3.0, maxiters=5, mask=xdf_image.mask)
        
            # find peaks
            sources_findpeaks = find_peaks(xdf_image.data, mask=xdf_image.mask, 
                                           threshold=10.*std, box_size=10,
                                           centroid_func=centroid_2dg)     
            # Display the table
            sources_findpeaks
            
   
            # plot the centroids of each of the sources
            
            # Set up the figure with subplots
            fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
            
            
            lower_value = np.percentile(data, lower_percentile)
            upper_value = np.percentile(data, upper_percentile)
            norm = Normalize(lower_value, upper_value)
            
            
            # Plot the data
            fitsplot = ax1.imshow(np.ma.masked_where(xdf_image.mask, xdf_image), norm=norm)
            ax1.scatter(sources_findpeaks['x_centroid'], sources_findpeaks['y_centroid'], s=10, marker='.', 
                        lw=1, alpha=0.7, color='r')#facecolor='None', edgecolor='r')
            
            # Define the colorbar
            cbar = plt.colorbar(fitsplot, fraction=0.046, pad=0.04, ticks=LogLocator(subs=range(10)))
            labels = ['$10^{-4}$'] + [''] * 8 + ['$10^{-3}$'] + [''] * 8 + ['$10^{-2}$']
            cbar.ax.set_yticklabels(labels)
            
            # Define labels
            cbar.set_label(r'Flux Count Rate ({})'.format(xdf_image.unit.to_string('latex')), 
                           rotation=270, labelpad=30)
            ax1.set_xlabel('X (pixels)')
            ax1.set_ylabel('Y (pixels)')
            ax1.set_title('find\_peaks Sources')

            
            
            # define the apertures
            circular_apertures = []
            for ap_row in range(len(sources_findpeaks)):
                position = (sources_findpeaks['x_centroid'][ap_row], sources_findpeaks['y_centroid'][ap_row])
                circular_apertures.append(CircularAperture(position, aperture_radius))
            
            """ eventually come up with a way to determine radii based on gaussians if it's not terribly difficult
                also, only make apertures around gaussian sources that aren't saturated"""
          
            
            
            
                        
            # Set up the figure with subplots
            fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
            
            # Plot the data
            fitsplot = ax1.imshow(np.ma.masked_where(xdf_image.mask, xdf_image), norm=norm)
            
            # Plot the apertures
            for aperture in circular_apertures:
                aperture.plot(color='red', alpha=0.7, axes = ax1)
            
            
            # Define the colorbar
            cbar = plt.colorbar(fitsplot, fraction=0.046, pad=0.04, ticks=LogLocator(subs=range(10)))
            labels = ['$10^{-4}$'] + [''] * 8 + ['$10^{-3}$'] + [''] * 8 + ['$10^{-2}$']
            cbar.ax.set_yticklabels(labels)
            
            # Define labels
            cbar.set_label(r'Flux Count Rate ({})'.format(xdf_image.unit.to_string('latex')), 
                           rotation=270, labelpad=30)
            ax1.set_xlabel('X (pixels)')
            ax1.set_ylabel('Y (pixels)')
            ax1.set_title('Circular Apertures')
            
            
            """ Performing Aperture Photometry """

            from photutils import aperture_photometry
            from astropy.table import QTable
            
            # use first aperture in image
            phot_datum = aperture_photometry(xdf_image, circular_apertures[0])
            phot_datum
                
            # The CCDData mask will be automatically applied
            phot_table = aperture_photometry(xdf_image, circular_apertures[0])
            id=1
            for aperture in circular_apertures[0:]:
                id += 1
                phot_row = aperture_photometry(xdf_image, aperture)[0]
                phot_row[0] = id
                phot_table.add_row(phot_row)
            
            # display the table
            phot_table 
            
           
           
          



#%%

