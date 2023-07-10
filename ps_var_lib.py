#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 12:12:02 2023

GOAL: 
    Maintain the functions for summer project

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
from astropy.nddata import CCDData
import imreg_dft

from photutils import DAOStarFinder

from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import ImageNormalize, LogStretch
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import numpy as np
from photutils.background import Background2D, MeanBackground





#%% make_image_user_input

def make_image_user_input(data, aperture, cmap = 'gray', plot_aperture = False,
               first_file = True, lower_percentile=0, upper_percentile=99):
    
    fig, ax = plt.subplots()
    
    
    scaled_data = np.sqrt(data)
    # scaled_data=data
    
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

def get_background_aperture(hdul, mask_sigma_coeff = 5):

    
    
    data = hdul[0].data
    header = hdul[0].header
    wcs = WCS(header)
    pixel_scale_y = header['CDELT2'] *u.deg
    
    # Convert the pixel scale to arcseconds
    pixel_scale = pixel_scale_y.to(u.arcmin).value


    # define the mask
    zero_mask = (data == 0)
    
    # store image 
    unit = u.electron/u.s
    image=CCDData(data, unit=unit, meta=header, mask=zero_mask)
    
    
    # stats
    median_counts = np.median(data) # calculate median counts per pixel of image
    stdev_counts = np.std(data)  # calculate standard deviation of the counts per pixel
    
    """ DETERMINE SUITABLE BACKGROUND REGION AND SUBTRACT FROM ENTIRE IMAGE """       
    
    radius=0.5 # radius units: arcmin
    
    aperture_radius = radius / pixel_scale # aperture radius has units of pixels
    
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
            
            # # apply the mask to the data (ignore pixels with counts > median+5stdev)
            # mask = np.where(data >= median_counts + mask_sigma_coeff*stdev_counts, True, False)
            # masked_data= np.ma.masked_array(data, mask)
            # # extract the pixel values within the aperture
            # mask = aperture.to_mask(method='center')
            # mask_image = mask.to_image(image.shape)
            # values = masked_data.data[np.where(mask_image)]
            
            # # calculate the median of unmasked pixels within the aperture
            # median_masked_counts = np.median(values)
    
            # Check if the current total counts are lower than the previous minimum
            if total_counts < min_median_counts:
                min_median_counts = total_counts
                best_center_x, best_center_y = center_x, center_y
    
    best_wcs_position = wcs.pixel_to_world(best_center_x, best_center_y) # convert to real world coords
    
    return best_wcs_position, aperture_radius



#%% process images

def process_images(fits_dir, mask_sigma_coeff, make_gif=False, cmap = 'magma'):
    """
    1) loop through fits files in directory
    2) read in important info (data, header etc)
    3) for the first file in directory, use get_background_aperture to find the 
        darkest circle of radius 0.5 arcsec - this region will be used to calculate
        background for all fits files in this directory
    4) change the value of infinite, nan, and negative pixel values to 0
    5) display image with aperture
    6) if make_gif = True, make a gif of the images plotted at the scale defined previously
    
    """
    
    
    # read in obj name and band
    split_path = fits_dir.split('/')
    name = split_path[6]
    band = split_path[8]
    
    # go through directory and get a path for each fits file 
    # don't work with the processed fits images (ending in _processed.fits), 
    # just the original files ending in mosaic-int.fits
    
    file_list = [file for file in os.listdir(fits_dir) if file.endswith('int.fits')]
    file_list_paths = []
    for file in file_list:
        file_list_paths.append(f'{fits_dir}/{file}')
        
        
    frames=[] # for gif
        
    # loop through and process files
    first_file = True # only define background aperture once, then change to False
    for image_file in file_list_paths:
        
        # open fits file
        hdul = fits.open(image_file)
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
        median_counts = np.median(data) # calculate median counts per pixel of image
        stdev_counts = np.std(data)  # calculate standard deviation of the counts per pixel
        
        # define the mask
        zero_mask = (data == 0)
        
        # store image 
        unit = u.electron/u.s
        image = CCDData(data, unit=unit, meta=header, mask=zero_mask)
        
        
        """ 
        Make background aperture at the darkest region of the first image, 
        and use this aperture to subtract the median (background) counts
        from the data. NOTE: mask is used to limit the impact of point sources 
        on our calculation of background brightness. 
        """
        if first_file == True:
            best_wcs_position, aperture_radius = get_background_aperture(hdul, mask_sigma_coeff)
        
        
        # use the median background counts in the background aperture as the background for the image
        # and subtract this background from all pixels
        
        # convert the ra and dec positions to pixels for the image   
        best_pixel_position = wcs.world_to_pixel_values(best_wcs_position.ra, best_wcs_position.dec)
        best_center_x = best_pixel_position[0]
        best_center_y = best_pixel_position[1]
            
        # Create the final aperture at the optimal center
        aperture = CircularAperture((best_center_x, best_center_y), aperture_radius)
        
        # apply the mask to the data (ignore pixels with counts > median+5stdev)
        mask = np.where(data >= median_counts + mask_sigma_coeff*stdev_counts, True, False)
        masked_data= np.ma.masked_array(data, mask)
        # extract the pixel values within the aperture
        mask = aperture.to_mask(method='center')
        mask_image = mask.to_image(image.shape)
        values = masked_data.data[np.where(mask_image)]
        
        # calculate the median of unmasked pixels within the aperture
        median_masked_counts = np.median(values)
        
        # Subtract the median background from the entire image
        data -= median_masked_counts
        
       
        """
        Make all infinite, negative, or nan pixel counts = 0 
        (do this after subtracting background so these pixels don't artificially
         impact our search for the darkest background aperture')
        """
        
        data[np.isnan(data)] = 0 
        data[np.isinf(data)] = 0
        data[data<0] = 0
       
        
        """
        Now, display an image with the background aperture plotted
        """
        
        if first_file == True:
            frame, lower_percentile, upper_percentile = \
                make_image_user_input(data=data, aperture=aperture,
                                      cmap=cmap, plot_aperture=True,
                                      first_file=True) 
            
            
        else: 
            frame, lower_percentile, upper_percentile = \
                make_image_user_input(data=data, aperture=aperture,
                                      cmap=cmap, plot_aperture=True,
                                      lower_percentile = lower_percentile, 
                                      upper_percentile = upper_percentile, 
                                      first_file = False) 

        
        
        frames.append(frame)
            
        # save the background_subtracted images as *_processed.fits 
        # new_file_name = image_file.replace('.fits', '_processed.fits')
        processed_fits = CCDData(data, unit=unit, meta=header)
        processed_fits.write(image_file.replace('.fits', '_processed.fits'), overwrite=True)    
            
        
        # hdul.writeto(os.path.join(fits_dir, new_file_name))
        hdul.close()
        
        first_file = False
    
    
    frames[0].save(f'/users/dzakaria/DATA/dzfiles/animations/{name}_{band}_{cmap}.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)
    

    


#%% point source photometry

def point_source_photometry(processed_image_dir, lower_percentile, upper_percentile, fwhm = 4, threshold=5):
    # list of fits files to go through
    file_list = [file for file in os.listdir(processed_image_dir) if file.endswith('processed.fits')]
    
    for file in file_list:
        
        path= f'{processed_image_dir}/{file}'
        fig, ax = plt.subplots()
        with fits.open(path) as hdul:
            data = hdul[0].data
            header = hdul[0].header 
        
        
        wcs = WCS(hdul[0].header)
        pixel_scale_y = header['CDELT2'] *u.deg
        
        # Convert the pixel scale to arcseconds
        # pixel_scale_x = pixel_scale_x.to(u.arcmin).value
        pixel_scale = pixel_scale_y.to(u.arcmin).value
        
        # # define the mask
        # mask = (data == 0)
        
        # store image 
        unit = u.electron/u.s
        xdf_image=CCDData(data, unit=unit, meta=header)
        
        

        """ CREATING COMP STAR APERTURES """
        
        from photutils import find_peaks
        from photutils.centroids import centroid_2dg
    
        # calculate statistics
        mean, median, std = sigma_clipped_stats(xdf_image.data, sigma=3.0, maxiters=5, mask=xdf_image.mask)

        # find peaks
        star_finder = DAOStarFinder(threshold, fwhm) 
        
        stars_found = star_finder(xdf_image.data)
        
        
        
   
        # plot the centroids of each of the sources
        
        # Set up the figure with subplots
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
        
        
        lower_value = np.percentile(data, lower_percentile)
        upper_value = np.percentile(data, upper_percentile)
        norm = Normalize(lower_value, upper_value)
        

        
        
        # Plot the stars found and apertures
        x = stars_found['xcentroid']
        y = stars_found['ycentroid']
        
        
        aperture_radius = 5
        
        positions = np.column_stack((x,y))
        apertures = CircularAperture(positions, r=aperture_radius)
        
   
        
        fitsplot = ax1.imshow(xdf_image, cmap='gray', norm=norm)
        ax1.scatter(x,y, marker='.', facecolor='r', edgecolor='r', s=20)


        # define the apertures
        circular_apertures = []
        for ap_row in range(len(stars_found)):
            position = (stars_found['xcentroid'][ap_row], stars_found['ycentroid'][ap_row])
            aperture = CircularAperture(position, aperture_radius)
            circular_apertures.append(aperture)
            aperture.plot(color='red', alpha=0.7, axes = ax1)
            
            
        for i, (xcoord, ycoord) in enumerate(zip(x,y)):
            ax1.text(xcoord, ycoord, str(i), color='red', fontsize=12, va='bottom', ha='left')
            
        plt.colorbar(fitsplot, ax=ax1)
        
        # display initial plot
        plt.show()
        
        
        # Ask the user if they are content with the apertures
        content = input("Are you content with these apertures? (yes/no): ")

        apertures_to_remove = []
        while content.lower() != "yes":
           
            # Ask the user which apertures to remove
            what_to_remove = input("Enter the indices of the apertures to remove (comma-separated): ")
            for i in what_to_remove.split(','):
                apertures_to_remove.append(int(i))
                
            # Re-plot the image without the removed apertures
            fig, ax2 = plt.subplots(1, 1, figsize=(8, 8))
            fitsplot = ax2.imshow(xdf_image, cmap='gray', norm=norm)
            
            
            # Remove the selected apertures
            circular_apertures = [aperture for idx, aperture in enumerate(circular_apertures) if idx not in apertures_to_remove]

            x = np.delete(x, apertures_to_remove)
            y = np.delete(y, apertures_to_remove)
            
            for aperture in circular_apertures:
                aperture.plot(color='red', alpha=0.7, axes=ax2)
            
            for i, (xcoord, ycoord) in enumerate(zip(x,y)):
                ax2.text(xcoord, ycoord, str(i), color='red', fontsize=12, va='bottom', ha='left') 
                
            plt.colorbar(fitsplot, ax=ax2)
            
            plt.show() 

            # Ask the user if they are content with the apertures
            content = input("Are you content with these apertures? (yes/no): ")
            
            
            if content.lower() == 'yes':
                break
 
        
        apertures_kept = []
        apertures_kept.extend([aperture for idx, aperture in enumerate(circular_apertures) if idx not in apertures_to_remove])
            
        return apertures_kept
    
    