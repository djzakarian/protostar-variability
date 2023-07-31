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
import math
import numpy as np
from astropy.io import fits
from photutils.aperture import ( aperture_photometry, CircularAperture, CircularAnnulus, 
                                EllipticalAperture, SkyCircularAperture, SkyCircularAnnulus, SkyEllipticalAperture)
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image
from astropy.wcs import WCS
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from astropy.nddata import CCDData
from photutils import DAOStarFinder, Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, Column, MaskedColumn
from  astropy.timeseries import TimeSeries
import re
from astropy.time import Time

    
#%% pixel_to_radec

def pixel_to_radec(x,y, wcs):
    """ converts pixel coordinates to RA and Dec for a wcs object
    x , y: (floats), coordinates
    wcs: (WCS), Astropy WCS object containing the coordinate transformation info
    
    returns: tuple with RA and Dec (in deg)
    """
    # create a skycoord object with the pixel coordinates and WCS object
    sky_coord = wcs.pixel_to_world(x,y)
    # extract ra and dec values from the object
    ra = sky_coord.ra.deg
    dec = sky_coord.dec.deg
    return ra, dec


#%% norm

def norm(data, lower_percentile=20, upper_percentile=99.5):
    data=np.sqrt(data)
    lower_value = np.percentile(data, lower_percentile)
    upper_value = np.percentile(data, upper_percentile)
    norm = Normalize(lower_value, upper_value)
    return norm()

#%% make_image_user_input

def make_image_user_input(data, apertures = [], cmap = 'turbo', plot_apertures = False,
               first_file = True, lower_value=0, upper_value=3, show_date=False, date='', 
               show_axes=False, add_scale=False,sqrt=True, scale_values=False):
    

    fig, ax = plt.subplots()
    
    if sqrt==True:
        scaled_data = np.sqrt(data)
    else:
        scaled_data = data

    
    fig.clear()
    # ax = fig.add_subplot(111)
    # ax.axis('off')
    
    

    
    if first_file == True:
        
        while True:
            try:
                
                if scale_values==False:
                    lower_percentile = float(input("Enter the lower percentile (0-100): "))
                    upper_percentile = float(input("Enter the upper percentile (0-100): "))
                    lower_value = np.percentile(scaled_data, lower_percentile)
                    upper_value = np.percentile(scaled_data, upper_percentile)
                    norm = Normalize(lower_value, upper_value)
        
                else:
                    lower_value = float(input("Enter the lower value: "))
                    upper_value = float(input("Enter the upper value: "))
                    norm = Normalize(lower_value, upper_value)
                    

                fig.clear()
                ax = fig.add_subplot(111)
                
                if show_axes==False:
                    ax.axis('off')
        
                im = ax.imshow(scaled_data, cmap=cmap, origin='lower', norm=norm)
        
                if plot_apertures==True:
                    for i, aperture in enumerate(apertures):
                        # aperture.plot(color='red', lw=1.5, alpha=1)
                        x, y = aperture.positions
                        ax.text(x, y, str(i), color='red', fontsize=9, 
                                ha='center', va='center', bbox=dict(facecolor='white', alpha=0.4))
                        # ax.text(x, y, str(i), color='red', fontsize=12, 
                        #         ha='center', va='center')
        
                # ax.set_xlim(0, data.shape[1])
                # ax.set_ylim(0, data.shape[0])
                
                ax.set_aspect('equal')
                
                # add the date in white font on the to left corner if show_date = True
                if show_date == True:
                    ax.text(0.02, 0.98, date, transform=ax.transAxes,
                            color='white', fontsize=12, ha='left', va='top', bbox=dict(facecolor='black', alpha=0.7))
                
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
            except Exception as e:
                print(e)
                pass
    else:
        
       
        # lower_value = np.percentile(scaled_data, lower_percentile)
        # upper_value = np.percentile(scaled_data, upper_percentile)
        norm = Normalize(lower_value, upper_value)

        fig.clear()
        ax = fig.add_subplot(111)
        # ax.axis('off')
        
        if show_axes==False:
            ax.axis('off')

        im = ax.imshow(scaled_data, cmap=cmap, origin='lower', norm=norm)

        if plot_apertures==True:
            for i, aperture in enumerate(apertures):
                # aperture.plot(color='red', lw=1.5, alpha=1)
                x, y = aperture.positions
                ax.text(x, y, str(i), color='red', fontsize=9, 
                        ha='center', va='center', bbox=dict(facecolor='white', alpha=0.4))
                # ax.text(x, y, str(i), color='red', fontsize=12, 
                #         ha='center', va='center')
            
        # ax.set_xlim(0, data.shape[1])
        # ax.set_ylim(0, data.shape[0])
        ax.set_aspect('equal')
        
        # add the date in white font on the to left corner if show_date = True
        if show_date == True:
            ax.text(0.02, 0.98, date, transform=ax.transAxes,
                    color='white', fontsize=12, ha='left', va='top', bbox=dict(facecolor='black', alpha=0.7))
            
            
        canvas = FigureCanvasTkAgg(fig, master=tk.Tk())
        canvas.draw()
        canvas.get_tk_widget().pack()

        image = Image.frombytes('RGB', canvas.get_width_height(), fig.canvas.tostring_rgb())

        plt.close(fig)

        plt.imshow(image)
        plt.show()
    
    return image, lower_value, upper_value




#%% get_background_aperture

def get_background_aperture(hdul):

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
          
            # Check if the current total counts are lower than the previous minimum
            if total_counts < min_median_counts:
                min_median_counts = total_counts
                best_center_x, best_center_y = center_x, center_y
    
    best_wcs_position = wcs.pixel_to_world(best_center_x, best_center_y) # convert to real world coords
    
    return best_wcs_position, aperture_radius

#%% background_subtraction (Background2D) (DONT THINK THIS WILL WORK)

def process_images_background2d(fits_file_path, mask, plot_backgroun=True, plot_background_subtracted_image=True):
    """ perform background subtraction for a fits file using Background2D.
        return the HDU List with background subtracted"""
    
# # read in obj name and band
# split_path = fits_dir.split('/')
# name = split_path[6]
# band = split_path[8]

# # go through directory and get a path for each fits file 
# # don't work with the processed fits images (ending in _processed.fits), 
# # just the original files ending in mosaic-int.fits

# file_list = [file for file in os.listdir(fits_dir) if file.endswith('int.fits')]
# file_list_paths = []
# for file in file_list:
#     file_list_paths.append(f'{fits_dir}/{file}')
        
        
#     # loop through and process files
#     for image_file in file_list_paths:
    
    

    
    # open fits file
    with fits.open(fits_file_path) as hdul:

        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
        
        
        if mask== 'mid_third':
            image_height, image_width = data.shape
            
            roi_start = image_width//5
            roi_end = 4* image_width//5
            
            mask = np.zeros_like(data)
            mask[roi_start:roi_end,:]=1
            
       
        # calculate stats of the data
        mean, median, std = sigma_clipped_stats(data)
        
        # create a background estimator 
        bkg_estimator = MedianBackground()
        
        # perform background estimation
        background = Background2D(data, (30, 30), filter_size=(2,2), 
                                  bkg_estimator=bkg_estimator, mask=mask, exclude_percentile=90)
        # subtract the background from the data
        background_subtracted_data = data - background.background
        # create a new HDU and save it as a processed fits file
        processed_hdul = hdul.copy()
        processed_hdul[0].data=background_subtracted_data
        
        # save fits file
        
    # 
    background_scaled = np.sqrt(background.background)
    sub_image_scaled = np.sqrt(background_subtracted_data)
        
    norm_val = norm(sub_image_scaled, lower_percentile=20, upper_percentile=99)
    # plot the background and background subtracted image
    fig, axs = plt.subplots(1,2,figsize=(10,5))
    axs[0].imshow(background_scaled, cmap='turbo', norm=norm_val)
    axs[0].set_title('Background')
    axs[1].imshow(sub_image_scaled, cmap='turbo', norm=norm_val)
    axs[1].set_title('Background Subtracted Image')
    plt.tight_layout()
    plt.show()
        
        
    return processed_hdul
    
        
        
#%% process images bgrnd ap

def process_images_bgrnd_ap(fits_dir, mask_sigma_coeff, make_gif=False, cmap = 'turbo'):
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
            best_wcs_position, aperture_radius = get_background_aperture(hdul)
        
        
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
        Now, display an image with the background aperture plotted if make_gif=True
        """
        if make_gif == True:
            if first_file == True:
                frame, lower_value, upper_value = \
                    make_image_user_input(data=data, aperture=aperture,
                                          cmap=cmap, plot_apertures=True,
                                          first_file=True) 
                
                
            else: 
                frame, lower_value, upper_value = \
                    make_image_user_input(data=data, aperture=aperture,
                                          cmap=cmap, plot_apertures=True,
                                          lower_value = lower_value, 
                                          upper_value = upper_value, 
                                          first_file = False) 
    
            
            
            frames.append(frame)
                
        # save the background_subtracted images as *_processed.fits 
        # new_file_name = image_file.replace('.fits', '_processed.fits')
        processed_fits = CCDData(data, unit=unit, meta=header)
        processed_fits.write(image_file.replace('.fits', '_processed.fits'), overwrite=True)    
            
        
        # hdul.writeto(os.path.join(fits_dir, new_file_name))
        hdul.close()
        
        first_file = False
    
    if make_gif == True:
        frames[0].save(f'/users/dzakaria/DATA/dzfiles/animations/{name}_{band}_{cmap}.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)
    
#%% make gif
def make_gif(fits_dir,cmap = 'turbo', epochs_tab_path='', endswith='_processed.fits', do_sqrt=False, scale_values=False):
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
    
    # open epochs_tab
    epochs_tab=QTable.read(epochs_tab_path)
    
    
    # go through directory and get a path for each fits file 
    # don't work with the processed fits images (ending in _processed.fits), 
    # just the original files ending in mosaic-int.fits
    
    file_list = [file for file in os.listdir(fits_dir) if file.endswith(endswith)]
    file_list_paths = []
    for file in file_list:
        file_list_paths.append(f'{fits_dir}/{file}')
        
    sorted_file_list_paths = sorted(file_list_paths)


    frames=[] # for gif
        
    # loop through and process files
    first_file = True # only define background aperture once, then change to False
    
    for image_file in sorted_file_list_paths:
        
        # open fits file
        hdul = fits.open(image_file)
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
        median_counts = np.median(data) # calculate median counts per pixel of image
        stdev_counts = np.std(data)  # calculate standard deviation of the counts per pixel
        
        # # define the mask
        # zero_mask = (data == 0)
        
        # store image 
        unit = u.electron/u.s
        image = CCDData(data, unit=unit, meta=header)
        
        # use epochs_tab to find the middle date of the observation
        
        # find the epoch in the filename
        pattern = 'size0_0667_ep(\d+)'  
        match = re.search(pattern, image_file)
        
        obj_epoch = math.floor(float(match.group(1)))
        
        
        # use name and epoch to find the epochs_tab entry 
        # this will tell us the date range, 
        # which can be used to find the middle of the observation window
        epochs_tab_row = epochs_tab[(epochs_tab['obj_name']==name) & (epochs_tab['obj_epoch']==obj_epoch)]
        mjd_obs1, mjd_obs2 = epochs_tab_row['mjd_obs1'], epochs_tab_row['mjd_obs2']
        
        
        
        # calculate the middle date of observation 
        obs_date = (mjd_obs1[0] + mjd_obs2[0]) / 2   
        
        # obs_date=obs_date.value
        
        
        obs_date = round(obs_date, 2)
        

        if first_file == True:
            frame, lower_value, upper_value = \
                make_image_user_input(data=data,
                                      cmap=cmap, plot_apertures=False,
                                      first_file=True, show_date=True, date=obs_date, 
                                      show_axes=False, add_scale=True, sqrt=do_sqrt, scale_values=scale_values) 
            
            
        else: 
            frame, lower_value, upper_value = \
                make_image_user_input(data=data, 
                                      cmap=cmap, plot_apertures=False,
                                      lower_value = lower_value, 
                                      upper_value = upper_value, 
                                      first_file = False, show_date=True, date=obs_date, 
                                      show_axes=False, add_scale=True, sqrt=do_sqrt, scale_values=scale_values) 

        
        
        frames.append(frame)
                
      
        
        # hdul.writeto(os.path.join(fits_dir, new_file_name))
        hdul.close()
        
        first_file = False
    
    endswith_base = endswith.replace('.fits', '')
    frames[0].save(f'/users/dzakaria/DATA/dzfiles/animations/{name}_{band}_{cmap}_{endswith_base}.gif', save_all = True, append_images=frames[1:], duration=500, loop=0)
    

    


#%% determine_point_source_apertures
def determine_point_source_apertures(file, processed_images_dir, lower_value = 0,
                                     upper_value = 3, fwhm = 5.5, threshold=5,
                                     aperture_radius = 8.25, inner_radius = 12, outer_radius= 18):
    
    """ all radii in arcsec """
    
    path= f'{processed_images_dir}/{file}'
    # fig, ax = plt.subplots()
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header 
      
    wcs = WCS(hdul[0].header)
    pixel_scale_y = header['CDELT2'] *u.deg
    
    # Convert the pixel scale to arcseconds
    # pixel_scale_x = pixel_scale_x.to(u.arcsec).value
    pixel_scale = pixel_scale_y.to(u.arcsec).value
    
    # store image 
    unit = u.electron/u.s
    xdf_image=CCDData(data, unit=unit, meta=header)
    

    """ CREATING COMP STAR APERTURES """

    # calculate statistics
    mean, median, std = sigma_clipped_stats(xdf_image.data, sigma=3.0, maxiters=5, mask=xdf_image.mask)

    # find peaks
    star_finder = DAOStarFinder(threshold, fwhm)
    stars_found = star_finder(xdf_image.data)
    
    # while True:
    #     if len(stars_found) > 25:
    #         fwhm += 1
    #         star_finder = DAOStarFinder(threshold, fwhm)
    #         stars_found = star_finder(xdf_image.data)
            
    #     elif len(stars_found) < 10:
    #         fwhm -= 1
    #         star_finder = DAOStarFinder(threshold, fwhm)
    #         stars_found = star_finder(xdf_image.data)
    #         break
    #     else:
    #         break

    # plot the centroids of each of the sources
    
    # Set up the figure with subplots
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
    
    # lower_value = np.percentile(data, lower_percentile)
    # upper_value = np.percentile(data, upper_percentile)
    norm = Normalize(lower_value, upper_value)
    
    # Plot the stars found and apertures
    x = stars_found['xcentroid']
    y = stars_found['ycentroid']      

    fitsplot = ax1.imshow(xdf_image, cmap='turbo', norm=norm)
    ax1.scatter(x,y, marker='.', facecolor='r', edgecolor='r', s=20)

    # convert the radii to pixel values
    aperture_radius = aperture_radius/pixel_scale
    inner_radius = inner_radius/pixel_scale
    outer_radius = outer_radius/pixel_scale

    # define the apertures
    circular_apertures = []
    circular_annuli = []
    for ap_row in range(len(stars_found)):
        position = (stars_found['xcentroid'][ap_row], stars_found['ycentroid'][ap_row])
        aperture = CircularAperture(position, aperture_radius)
        circular_apertures.append(aperture)
        aperture.plot(color='red', alpha=0.7, axes = ax1)
        
        annulus_aperture = CircularAnnulus(position, r_in=inner_radius, r_out=outer_radius)
        circular_annuli.append(annulus_aperture)
        annulus_aperture.plot(color='blue', alpha=0.7, axes=ax1)
        
        
    for i, (xcoord, ycoord) in enumerate(zip(x,y)):
        ax1.text(xcoord, ycoord, str(i), color='red', fontsize=12, va='bottom', ha='left')
        
    plt.colorbar(fitsplot, ax=ax1)
    
    plt.show() # display initial plot
    
    # # Ask the user if they are content with the apertures
    # content = input("Are you content with these apertures? (yes/no): ")

    # apertures_to_remove = []
    # while content.lower() != "yes":
       
    #     # Ask the user which apertures to remove
    #     what_to_remove = input("Enter the indices of the apertures to remove (comma-separated): ")
    #     for i in what_to_remove.split(','):
    #         apertures_to_remove.append(int(i))
            
    #     # Re-plot the image without the removed apertures
    #     fig, ax2 = plt.subplots(1, 1, figsize=(8, 8))
    #     fitsplot = ax2.imshow(xdf_image, cmap='gray', norm=norm)
        
        
    #     # Remove the selected apertures
    #     circular_apertures = [aperture for idx, aperture in enumerate(circular_apertures) if idx not in apertures_to_remove]

    #     x = np.delete(x, apertures_to_remove)
    #     y = np.delete(y, apertures_to_remove)
        
    #     for aperture in circular_apertures:
    #         aperture.plot(color='red', alpha=0.7, axes=ax2)
        
    #     for i, (xcoord, ycoord) in enumerate(zip(x,y)):
    #         ax2.text(xcoord, ycoord, str(i), color='red', fontsize=12, va='bottom', ha='left') 
            
    #     plt.colorbar(fitsplot, ax=ax2)
        
    #     plt.show() 

    #     # Ask the user if they are content with the apertures
    #     content = input("Are you content with these apertures? (yes/no): ")
        
    #     if content.lower() == 'yes':
    #         break

    # apertures_kept = []
    # apertures_kept.extend([aperture for idx, aperture in enumerate(circular_apertures) if idx not in apertures_to_remove])
        

    # convert the apertures to RA and Dec
    wcs_apertures = []
    wcs_annuli=[]
    for aperture in circular_apertures:

        xcenter, ycenter = aperture.positions
        ra, dec = pixel_to_radec(xcenter, ycenter, wcs)
        center_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
        wcs_aperture =  SkyCircularAperture(center_coord, aperture_radius * u.pixel)
        wcs_apertures.append(wcs_aperture)
        
        wcs_annulus = SkyCircularAnnulus(center_coord, r_in=inner_radius*u.pixel,
                                          r_out=outer_radius*u.pixel)
        wcs_annuli.append(wcs_annulus)
        
    
    return wcs_apertures, wcs_annuli

#%% create_elliptical_apertures

def create_elliptical_aperture(x, y, semi_major_axis, semi_minor_axis, theta):
    position = (x, y)
    a = semi_major_axis
    b = semi_minor_axis
    aperture = EllipticalAperture(position, a, b, theta)
    return aperture



#%% determine target apertures

def determine_target_apertures(combined_image_file, processed_images_dir, fwhm = 5.5, threshold=5,
                                     aperture_radius = 8.25, inner_radius = 12, outer_radius= 18):
    
    path= f'{processed_images_dir}/{combined_image_file}'
    # fig, ax = plt.subplots()
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header 
      
    wcs = WCS(hdul[0].header)
    pixel_scale_y = header['CDELT2'] *u.deg
    
    # Convert the pixel scale to arcseconds
    # pixel_scale_x = pixel_scale_x.to(u.arcsec).value
    pixel_scale = pixel_scale_y.to(u.arcsec).value
    
    # store image 
    unit = u.electron/u.s
    xdf_image=CCDData(data, unit=unit, meta=header)
    

    """ CREATING APERTURE(S) """
    
    image, lower_value, upper_value = make_image_user_input(data)
    
    apertures = []
    while True:
                
        try:
            print(f"\nDefining Aperture:")
    
            x = float(input("aperture x coord: "))
            y = float(input("aperture y coord: "))
            semi_major_axis = float(input("semi-major axis: "))
            semi_minor_axis = float(input("semi-minor axis: "))
            theta = float(input("angle (deg): "))
            theta = theta * math.pi/180
            aperture = create_elliptical_aperture(x, y, semi_major_axis, semi_minor_axis, theta)
            
            test_aperture =[aperture]
            for ap in apertures:
                test_aperture.append(ap)
    
            # plot the aperture to visualize it
            image, lower, upper = make_image_user_input(data, apertures = test_aperture, 
                                  plot_apertures=True, upper_value=upper_value,
                                  lower_value = lower_value, first_file=False)
            
            response = input("\nAre you happy with this aperture? (yes/no): ").lower()
            if response == 'yes':
                apertures.append(aperture)
            
        
            response = input("\nAre you done adding apertures?: ").lower()
            if response == 'yes': 
                
                exit_loop = True
                # image file name
                image_file = combined_image_file.replace('coadd', 'target_ap')
                image_file = image_file.replace('.fits', '.jpg')
                image.save(f'{processed_images_dir}/{image_file}')
                break
        
        
        except Exception as e:
            print(e)
            continue         
        
    
    

    # convert the apertures to RA and Dec
    target_wcs_apertures = []
    for aperture in apertures:
        a = aperture.a
        b = aperture.b
        theta = aperture.theta
        xcenter, ycenter = aperture.positions
        ra, dec = pixel_to_radec(xcenter, ycenter, wcs)
        center_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
        target_wcs_aperture =  SkyEllipticalAperture(center_coord, a * u.pixel, b*u.pixel, theta*u.rad)
        target_wcs_apertures.append(target_wcs_aperture)
        
       

    
    return target_wcs_apertures
   
#%% process_images_old

def process_images(fits_dir, mask_sigma_coeff, make_gif=False, cmap = 'turbo'):
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
        
#%% time series data old

def time_series_point_source_photometry(processed_image_dir):
    # list of fits files to go through
    file_list = [file for file in os.listdir(processed_image_dir) if file.endswith('processed.fits')]
    
    # first determine the apertures from the first image in the directory
    first_file = file_list[0]
    apertures = determine_point_source_apertures(first_file, processed_image_dir)
    
    

    
    # Create an empty table to store the photometry results
    phot_table = Table(names=['File', 'Date', 'Aperture', 'Flux'], dtype=['U50', 'U50', int, float])
    
    # Iterate over each FITS file
    for file in file_list:
        # Read the FITS file
        path = os.path.join(processed_image_dir, file)
        with fits.open(path) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            date = header.get('DATE-OBS', 'Unknown')
    
        # Perform aperture photometry using the same set of apertures
        positions = [(aperture.ra.deg, aperture.dec.deg) for aperture in apertures]
        aperture_radius = 5  # Adjust as needed
        apertures = CircularAperture(positions, r=aperture_radius)
        phot_table_per_file = aperture_photometry(data, apertures)
    
        # Extract the flux values
        fluxes = phot_table_per_file['aperture_sum']
    
        # Add the results to the main photometry table
        for i, flux in enumerate(fluxes):
            phot_table.add_row([file, date, i, flux])
    
    # Save the photometry table to a file
    output_file = '{processed_image_dir}/photometry_table.fits'
    phot_table.write(output_file, overwrite=True)
    
    # Display the photometry table
    print(phot_table)



#%% time series tab

def time_series_tab(processed_images_dir, epochs_tab_path, 
                    wcs_apertures, wcs_annuli, target_wcs_apertures):
    
    """
    make a time series table for a given directory of fits files of a target
    
    args:
        processed_images_dir: path where the fits files are stored
        epochs_tab_path: path to epochs_table with correct epochs infomration
        wcs_apertures: list of photutils wcs aperture objects
    """
    
    # open epochs_tab
    epochs_tab=QTable.read(epochs_tab_path)
    
    time_list = []
    table = Table()
    
    # table_data = {'Time': []}
    table.meta['description'] = 'this table provides the calibrated magnitudes for each aperture'
    
    for row in range(len(target_wcs_apertures)):
        
        
        aperture = target_wcs_apertures[row]
        

        ra = aperture.positions.ra.deg
        dec = aperture.positions.dec.deg
        ap_a = aperture.a
        ap_b = aperture.b
        ap_theta = aperture.theta
        
        new_column  = MaskedColumn(name=f'TARGET aperture:{row}', dtype=float)
        table.add_column(new_column)
        
        table[f'TARGET aperture:{row}'].meta['ra[deg]'] = ra
        table[f'TARGET aperture:{row}'].meta['dec[deg]'] = dec
        table[f'TARGET aperture:{row}'].meta['semi-major_axis[pix]'] = ap_a
        table[f'TARGET aperture:{row}'].meta['semi-minor_axis[pix]'] = ap_b
        table[f'TARGET aperture:{row}'].meta['theta[rad]'] = ap_theta
        
         
    for row in range(len(wcs_apertures)):
        aperture = wcs_apertures[row]
        annulus = wcs_annuli[row]
        ra = aperture.positions.ra.deg
        dec = aperture.positions.dec.deg
        ap_rad = aperture.r
        an_inner_r = annulus.r_in
        an_outer_r = annulus.r_out
        
        new_column  = MaskedColumn(name=f'aperture:{row + len(target_wcs_apertures)}', dtype=float)
        table.add_column(new_column)
        
        table[f'aperture:{row + len(target_wcs_apertures)}'].meta['ra[deg]'] = ra
        table[f'aperture:{row + len(target_wcs_apertures)}'].meta['dec[deg]'] = dec
        table[f'aperture:{row + len(target_wcs_apertures)}'].meta['aperture_rad[pix]'] = ap_rad
        table[f'aperture:{row + len(target_wcs_apertures)}'].meta['annulus_inner_rad[pix]'] = an_inner_r
        table[f'aperture:{row + len(target_wcs_apertures)}'].meta['annulus_outer_rad[pix]'] = an_outer_r
     
        
       
    
    file_list = [file for file in os.listdir(processed_images_dir) if file.endswith('_processed.fits')]
    file_list_paths = []
    for file in file_list:
        file_list_paths.append(f'{processed_images_dir}/{file}')
        
    sorted_file_list_paths = sorted(file_list_paths)
    
    time_row_counter = 0
    
    for file in sorted_file_list_paths:
        
        """ 
        file: original, processed image file path
        sub_file: associated image path after coadd subtraction
        div_file: associated image path after coadd subtraction
        """
        sub_file = file.replace('_processed.fits', '_processed_subtracted.fits')
        div_file = file.replace('_processed.fits', '_processed_divided.fits')
        
        table.add_row()
        
        # find the epoch in the filename
        pattern = 'size0_0667_ep(\d+)'
        
        
        match = re.search(pattern, file)
        
        obj_epoch = math.floor(float(match.group(1)))
        
        # find obj_name 
        path_components = file.split('/')
        obj_name = path_components[6]
        
        # use name and epoch to find the epochs_tab entry 
        # this will tell us the date range, 
        # which can be used to find the middle of the observation window
        epochs_tab_row = epochs_tab[(epochs_tab['obj_name']==obj_name) & (epochs_tab['obj_epoch']==obj_epoch)]
        mjd_obs1, mjd_obs2 = epochs_tab_row['mjd_obs1'], epochs_tab_row['mjd_obs2']
        
        # calculate the middle date of observation 
        time = (mjd_obs1 + mjd_obs2) / 2
        
        time_list.append(time)
    
        with fits.open(file) as hdul:
            wcs = WCS(hdul[0].header)
            data = hdul[0].data
            
            # loop through the targets
            for row in range(len(target_wcs_apertures)):
                header = hdul[0].header
                sky_aperture = target_wcs_apertures[row]
                aperture_pixels = sky_aperture.to_pixel(wcs)
                phot_table = aperture_photometry(hdul[0].data, [sky_aperture], wcs=wcs)
                flux = phot_table['aperture_sum_0']
                

                # now convert to magnitudes 
                # (https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#WISEPS)
                
                # magnitude zero point from header
                magzp = header['MAGZP']
                
                # aperture correction given by WISE
                # (https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4c.html#circ)
                band = header['BAND']
                if band == 1:
                    ac = 0.222
                elif band == 2:
                    ac = 0.280
                
                # if flux is negative, ignore the log term
                if flux >= 0:
                    mag = magzp - 2.5 * np.log10(flux)
                
                else:
                    mag = magzp
                    
                    
                table[f'TARGET aperture:{row}'][time_row_counter] = mag
                
            
            # loop through the comp stars
            for row in range(len(wcs_apertures)):
                header = hdul[0].header
                sky_aperture = wcs_apertures[row]
                sky_annulus = wcs_annuli[row]
                aperture_pixels = sky_aperture.to_pixel(wcs)
                annulus_pixels = sky_annulus.to_pixel(wcs)
                phot_table = aperture_photometry(data, [sky_aperture, sky_annulus], wcs=wcs)
                aperture_flux = phot_table['aperture_sum_0']
                annulus_flux = phot_table['aperture_sum_1']
                bg_flux = annulus_flux * aperture_pixels.area/annulus_pixels.area
                flux = aperture_flux - bg_flux
                
                # now convert to magnitudes 
                # (https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#WISEPS)
                
                # magnitude zero point from header
                magzp = header['MAGZP']
                
                # aperture correction given by WISE
                # (https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4c.html#circ)
                band = header['BAND']
                if band == 1:
                    ac = 0.222
                elif band == 2:
                    ac = 0.280
                
                # if flux is negative, ignore the log term
                if flux >= 0:
                    mag = magzp - 2.5 * np.log10(flux) - ac
                
                else:
                    mag = magzp - ac
                    
                    
                table[f'aperture:{row + len(target_wcs_apertures)}'][time_row_counter] = mag
                
                
                # pixel_aperture = CircularAperture(aperture_pixels.positions, aperture_pixels.r)
                # aperture_sum = aperture_photometry(hdul[0].data, pixel_aperture)['aperture_sum']

                # table[f'aperture:{row}']
                #     table_data[f'({aperture.positions.ra.deg}, {aperture.positions.dec.deg})'].append(aperture_sum)
                 

        time_row_counter += 1
        
    time = np.concatenate(time_list)
    
    

    # for aperture in wcs_apertures:
    #     # flux = np.concatenate(flux_dict[aperture])
    
    time = Time(time, format='mjd')
    
    return TimeSeries(table, time=time)

            
            
#%% plot lightcurve

def plot_lightcurve(time_series_tab, obj_name, band,endswith, 
                    style='seaborn-muted', y_lower=0, y_upper=15, 
                    just_target=False, show_plot=False):
    
    plt.style.use(style)
    x=np.array(time_series_tab['time'].value)
    labels=[]
    
    band=band.replace('b','W')
    # determine title of the plot
    if just_target==False:
        title_name = f'{obj_name} {band}'
    elif endswith == '_processed.fits':
        title_name = f'{obj_name} {band}' 
        
        
    for column in time_series_tab.colnames:
        
        if just_target == True:
            
            if column != 'time' and 'TARGET' in column:
                y = np.array(time_series_tab[column])
                plt.scatter(x,y, marker='o',s = 15)
                plt.ylim(y_lower, y_upper)
                labels.append(column)
                plt.xlabel('mjd date')
                plt.ylabel('magnitude')
                plt.title(title_name)
        else: # plot compstars   

            if column != 'time':
                
                # if the metadata has been updated
                try:
                    if 'TARGET' in column:
                        y = np.array(time_series_tab[column])
                        plt.scatter(x,y, marker='o', s = 15)
                        plt.ylim(y_lower, y_upper)
                        labels.append(column)
                    elif time_series_tab[column].meta['ignore_aperture']=='False':
                        y = np.array(time_series_tab[column])
                        plt.scatter(x,y, marker='D', s = 6)
                        plt.ylim(y_lower, y_upper)
                        labels.append(column)
                        plt.xlabel('mjd date')
                        plt.ylabel('magnitude')
                        plt.title(title_name)

                        
                except Exception as e:
                    # # checkpoint
                    # print(f'exception1:{e}')
                    y = np.array(time_series_tab[column])
                    plt.scatter(x,y, s = 6)
                    plt.ylim(y_lower, y_upper)
                    labels.append(column)
                    plt.xlabel('mjd date')
                    plt.ylabel('magnitude')
                    plt.title(title_name)
                    
    plt.legend(labels, loc='upper left', ncol=3, fontsize='small')
    
    if show_plot==True:
        plt.show()
    else:
        return plt
    
   
    
    
    

#%% make light curve
def make_light_curve(fits_dir,cmap = 'gray', endswith='_processed.fits', 
                     do_sqrt=False, scale_values=False, just_target=False, 
                     style='seaborn-muted', lower_value=0, upper_value=3):
    """
    1) open the time series table
    2) display the coadd image with apertures plotted
    3) choose which apertures to keep
    4) choose y min and maxes for plots
    5) save plot
    
    """
    plt.style.use(style)
    
    # read in obj name and band
    split_path = fits_dir.split('/')
    name = split_path[6]
    band = split_path[8]
    
    time_series_tab_path = f'{fits_dir}/time_series_{name}_{band}.ecsv'
    
    endswith_base = endswith.replace('.fits', '')
    # open time_series_tab
    time_series_tab=QTable.read(time_series_tab_path)
    
    # open the coadd image 
    image_file = f'{fits_dir}/coadd_{name}_{band}.fits'
    coadd_hdul = fits.open(image_file)
    coadd_data = coadd_hdul[0].data
    coadd_header = coadd_hdul[0].header
    coadd_wcs = WCS(coadd_header)

    
    colnames = time_series_tab.colnames
    colnames.remove('time')

    
    apertures=[]
    for ap_name in colnames:
        if "TARGET" in ap_name:
            ra = time_series_tab[ap_name].meta['ra[deg]']
            dec = time_series_tab[ap_name].meta['dec[deg]']
            a = time_series_tab[ap_name].meta['semi-major_axis[pix]'].value
            b = time_series_tab[ap_name].meta['semi-minor_axis[pix]'].value
            theta = time_series_tab[ap_name].meta['theta[rad]'].value
            
            aperture_position = coadd_wcs.world_to_pixel_values(ra, dec)
            # x = aperture_position[0]
            # y = aperture_position[1]

            aperture = EllipticalAperture(aperture_position, a, b, theta)
            apertures.append(aperture)
        else: 
            ra = time_series_tab[ap_name].meta['ra[deg]']
            dec = time_series_tab[ap_name].meta['dec[deg]']
            rad = time_series_tab[ap_name].meta['aperture_rad[pix]'].value
            aperture_position = coadd_wcs.world_to_pixel_values(ra, dec)
            aperture = CircularAperture(aperture_position, rad)
            apertures.append(aperture)
    
    l, u = lower_value, upper_value
    if just_target==False:
        image, l, u = make_image_user_input(coadd_data, apertures=apertures, plot_apertures=True, cmap=cmap)
    else:
        make_image_user_input(coadd_data, first_file=False,
                              apertures=apertures, plot_apertures=True, 
                              cmap=cmap, lower_value=l, upper_value=u)
    
    while True:
        try:
                
            apertures_to_keep = []
    
    
            y_lower = float(input("Enter the lower y limit: "))
            y_upper = float(input("Enter the upper y limit: "))
            plot_lightcurve(time_series_tab, obj_name=name, band=band,endswith=endswith, 
                            y_lower=y_lower, y_upper=y_upper, style=style, show_plot=True)
                    
            if just_target == False:
                # Ask the user which apertures to remove
                what_to_keep = input("Enter the indices of the apertures to keep (comma-separated): ")
                for i in what_to_keep.split(','):
                    apertures_to_keep.append(int(i))
                
                
                for ap_name in colnames:
                    ap = ap_name.replace('TARGET aperture:','')
                    ap = ap.replace('aperture:','')
                    # print(f'ap:{ap}, aps to keep:{apertures_to_keep}')
                    if int(ap) in apertures_to_keep:
                        time_series_tab[ap_name].meta['ignore_aperture']= 'False'
                        
                        # checkpoint
                        # print('checkpoint',time_series_tab[ap_name].meta['ignore_aperture'])
                    else:
                        time_series_tab[ap_name].meta['ignore_aperture'] = 'True'
    

            plot_lightcurve(time_series_tab,obj_name=name, band=band, 
                            endswith=endswith, y_lower=y_lower, y_upper=y_upper, 
                            just_target=just_target, show_plot=True)
            
            kept_apertures=[]
            for i, aperture in enumerate(apertures):
                if i in apertures_to_keep:
                    kept_apertures.append(aperture)
                
            image, l, u = make_image_user_input(coadd_data, apertures=kept_apertures, 
                                  first_file=False,lower_value=l, upper_value=u, plot_apertures=True, cmap=cmap)
                    
                    
            content = input("Are you content with the plot? (yes/no): ")
            
            
            if content.lower() == "yes":
                if just_target==True:
                    
                   
                    # remake final image to save it
                    plot = plot_lightcurve(time_series_tab,obj_name=name, band=band, 
                                          endswith=endswith, y_lower=y_lower, y_upper=y_upper,
                                          just_target=just_target, show_plot=False)
                    
                    plot_save_path = f'/users/dzakaria/DATA/dzfiles/lightcurves/{name}_{band}{endswith_base}_target_lc.png'
                    plot.savefig(plot_save_path)
                    
                    
                else: # if compstars are also plotted
                
                    # remake final image to save it
                    plot = plot_lightcurve(time_series_tab,obj_name=name, band=band, 
                                          endswith=endswith, y_lower=y_lower, y_upper=y_upper,
                                          just_target=just_target, show_plot=False)
                    
                
                    plot_save_path = f'/users/dzakaria/DATA/dzfiles/lightcurves/{name}_{band}{endswith_base}_target_and_compstars_lc.png'
                    ap_save_path = f'/users/dzakaria/DATA/dzfiles/lightcurves/{name}_{band}{endswith_base}_target_and_compstars_aps.png'
                    plt.savefig(plot_save_path)
                    image.save(ap_save_path, format="PNG")
                break
        except Exception as e:
            # checkpoint
            print(f'exception2:{e}')
            pass

    
    coadd_hdul.close()
    time_series_tab.write(time_series_tab_path, format='ascii.ecsv', overwrite=True)
    
    return l, u


#%% combine_images

def combine_images(processed_images_dir):
    file_list = [file for file in os.listdir(processed_images_dir) if file.endswith('_processed.fits')]
    file_list_paths = []
    for file in file_list:
        file_list_paths.append(f'{processed_images_dir}/{file}')
        
    
    # open the first fits file to get header info
    # had to change it from the first to the last because 
    # the file isn't closing properly before looping through the directory
    first_file = fits.open(file_list_paths[0])
    header = first_file[0].header
    
    # print(header)
    
    # get the shape of the data arrays
    data_shape = first_file[0].data.shape
    
    # create an empty array to store the co-added data
    coadd_data = np.zeros(data_shape)
    
    first_file.close() # close the first fits file 
    

    # loop over the fits files and add their data arrays to the coadd array
    for fits_file in file_list_paths:
        
        hdulist = fits.open(fits_file)
        data = hdulist[0].data
        coadd_data += data
        hdulist.close()
        
    # create a new fits hdu object with the coadded data and header
    coadd_hdu = fits.PrimaryHDU(data=coadd_data, header=header)
    
    # save the coadded file
    split_path = processed_images_dir.split('/')
    name = split_path[6]
    band = split_path[8]
    coadd_hdu.writeto(f'{processed_images_dir}/coadd_{name}_{band}.fits', overwrite=True)
        
    # make_image_user_input(coadd_data) # view a plot of the image
    
    
    
#%% mask_target
def mask_target(processed_images_dir, sigma_coeff=1):
    # open file
    image_file = f'{processed_images_dir}/coadd.fits'
    
    hdul = fits.open(image_file)
    data = hdul[0].data
    header = hdul[0].header
    wcs = WCS(header)
    median_counts = np.median(data) # calculate median counts per pixel of image
    stdev_counts = np.std(data)  # calculate standard deviation of the counts per pixel
    # pixel_scale_y = header['CDELT2'] *u.deg
    # # Convert the pixel scale to arcseconds
    # pixel_scale = pixel_scale_y.to(u.arcmin).value
    
    # set all pixels with counts less than the median + 3 sigma = 0 
    data[data < median_counts + sigma_coeff* stdev_counts] = 0
    
    make_image_user_input(data)
    
    
    hdul.close()
    
#%% subtract_coadd
def subtract_coadd(processed_images_dir):
    
    # save the coadded file
    split_path = processed_images_dir.split('/')
    name = split_path[6]
    band = split_path[8]

    image_file = f'{processed_images_dir}/coadd_{name}_{band}.fits'
    
    coadd_hdul = fits.open(image_file)
    coadd_data = coadd_hdul[0].data
    coadd_header = coadd_hdul[0].header
    coadd_wcs = WCS(coadd_header)
    median_counts = np.median(coadd_data) # calculate median counts per pixel of image
    stdev_counts = np.std(coadd_data)  # calculate standard deviation of the counts per pixel
    
    # print(median_counts)

    
    """ loop through the processed files and subtract the avg from the epoch image """
    
    # count the number of files
    file_counter = 0
    
    file_list = [file for file in os.listdir(processed_images_dir) if file.endswith('processed.fits')]
    file_list_paths = []
    for file in file_list:
        file_list_paths.append(f'{processed_images_dir}/{file}')
        file_counter += 1
            
   
    for image_file in file_list_paths:
        
        # open fits file
        hdul = fits.open(image_file)
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
        median_counts = np.median(data) # calculate median counts per pixel of image
        stdev_counts = np.std(data)  # calculate standard deviation of the counts per pixel

        scale_ratio = 1 / (file_counter)
        
       
        sub_data = data - coadd_data * ( scale_ratio )
        
        # print(f'coadd_data: {coadd_data}')
        # print(f'\ndata: {data}')
        # print(f'\nsub_data: {sub_data}')
        
        # create a new fits hdu object with the coadded data and header
        sub_hdu = fits.PrimaryHDU(data=sub_data, header=header)
        
        # new file name 
        sub_file_name =  image_file.replace('_processed.fits', '_processed_subtracted.fits')
        
        # save the coadded file
        sub_hdu.writeto(f'{sub_file_name}', overwrite=True)
            
        # make_image_user_input(sub_data) # view a plot of the image 
    
    
        hdul.close()
        coadd_hdul.close()
        

#%% divid_coadd
def divide_coadd(processed_images_dir):
    
    # save the coadded file
    split_path = processed_images_dir.split('/')
    name = split_path[6]
    band = split_path[8]

    image_file = f'{processed_images_dir}/coadd_{name}_{band}.fits'
    
    coadd_hdul = fits.open(image_file)
    coadd_data = coadd_hdul[0].data
    coadd_header = coadd_hdul[0].header
    coadd_wcs = WCS(coadd_header)
    median_counts = np.median(coadd_data) # calculate median counts per pixel of image
    stdev_counts = np.std(coadd_data)  # calculate standard deviation of the counts per pixel
    
    # print(median_counts)

    
    """ loop through the processed files and subtract the avg from the epoch image """
    
    # count the number of files
    file_counter = 0
    
    file_list = [file for file in os.listdir(processed_images_dir) if file.endswith('processed.fits')]
    file_list_paths = []
    for file in file_list:
        file_list_paths.append(f'{processed_images_dir}/{file}')
        file_counter += 1
            
   
    for image_file in file_list_paths:
        
        # open fits file
        hdul = fits.open(image_file)
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
        median_counts = np.median(data) # calculate median counts per pixel of image
        stdev_counts = np.std(data)  # calculate standard deviation of the counts per pixel

        scale_ratio = 1 / (file_counter)
        
       
        sub_data = data / (coadd_data * ( scale_ratio ))
        
        # print(f'coadd_data: {coadd_data}')
        # print(f'\ndata: {data}')
        # print(f'\nsub_data: {sub_data}')
        
        # create a new fits hdu object with the coadded data and header
        sub_hdu = fits.PrimaryHDU(data=sub_data, header=header)
        
        # new file name 
        sub_file_name =  image_file.replace('_processed.fits', '_processed_divided.fits')
        
        # save the coadded file
        sub_hdu.writeto(f'{sub_file_name}', overwrite=True)
            
        # make_image_user_input(sub_data) # view a plot of the image 
    
    
        hdul.close()
        coadd_hdul.close()
        
#%% pipeline
      
def pipeline(file_dir, epochs_tab_path , process_files = True, comp_star_phot = True, combine = True,
             subtract = True, divide = True, point_source_fwhm = 5.5, point_source_threshold = 5, 
             compstar_ylim = 5000, target_wcs_apertures=[], wcs_apertures=[], wcs_annuli = []):
    
    
    """ file directory is the directory where the fits files for a given object + pipeline are stored.
        epochs_tab_path is the path where the epochs table information is stored.
        process_files, comp_star_phot, combine, subtract, and divide all correspond to parts of the pipeline.
        Set these to 'False' if they have already been done for a target."""
    

    # get info from file directory
    split_path = file_dir.split('/')
    name = split_path[6]
    band = split_path[8]
    
    
   
    if process_files == True:
        """ subtract the background from the images (if they haven't already been processed) """
        # and save the new fits files with the suffix _processed.fits 
        process_images_bgrnd_ap(file_dir, mask_sigma_coeff = 5, cmap = 'turbo', make_gif=False)
      
   
    if comp_star_phot == True:
        
        """ use the first fits file to determine apertures for the comp stars
            it may choose some bad apertures but they can be ignored.
            also, the apertures are found for b1 and used for both bands"""   
        
        
        if band == 'b1':
            file_list = [file for file in os.listdir(file_dir) if file.endswith('_processed.fits') and 'ep0' in file]
            file_list_paths = []
            for file in file_list:
                file_list_paths.append(f'{file_dir}/{file}')
              
            # combined image path 
            combined_file_path = f'coadd_{name}_{band}.fits'
            
            target_wcs_apertures = determine_target_apertures(combined_file_path, file_dir)
            wcs_apertures, wcs_annuli = determine_point_source_apertures(file_list[0], file_dir, fwhm = point_source_fwhm, threshold = point_source_threshold)
            # if band = b2, then the apertures and annuli are read into the function after being previously found
        
        """ make a time series table using all of the apertures found previously - 
            the aperture and annulus info are saved in the metadata of the table """
          
        
        ts_tab= time_series_tab(file_dir, epochs_tab_path, wcs_apertures, wcs_annuli, target_wcs_apertures)
        # update metadata for the table with the fwhm and threshold used
        ts_tab.meta['point_source_selection_parameters'] = f'fwhm={point_source_fwhm}, threshold={point_source_threshold}'
        # save the file in the directory
        ts_tab.write(f'{file_dir}/time_series_{name}_{band}.ecsv', format='ascii.ecsv', overwrite=True)
        ts_tab.write(f'{file_dir}/time_series_{name}_{band}.csv', format='ascii.csv', overwrite=True)
        
        # """ plot the lightcurves of the compstars """
        # plot_lightcurve(ts_tab, y_upper = compstar_ylim)
    
    if combine == True:
        """ combine all of the images into an average image """
        combine_images(file_dir)
    
    if subtract == True:
        """ scale and subtract coadd from each image - 
            each image will be saved with the suffix _subtracted
            (so complete suffix will be _processed_subtracted.fits )"""
            
        subtract_coadd(file_dir)
    
    if divide == True:
        """ scale and divide coadd from each image -
            each image will be saved with the suffix _divided """
        divide_coadd(file_dir)
    
    return wcs_apertures, wcs_annuli, target_wcs_apertures
    
#%% run_pipeline
def run_pipeline(coadd_directory, obj_names, bands, epochs_tab_path) :
    processed_names = []

    for name in obj_names:
        if name in processed_names:
            continue
        for band in bands:
            print(name)
            
            file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
            
            if band == 'b1':
                wcs_apertures, wcs_annuli, target_wcs_apertures = pipeline(file_dir, epochs_tab_path = epochs_tab_path)
                
            else:
                pipeline(file_dir, epochs_tab_path = epochs_tab_path, wcs_apertures = wcs_apertures, 
                         wcs_annuli = wcs_annuli, target_wcs_apertures = target_wcs_apertures)
                
                # after successfully completed, add to list of processed names
                processed_names.append(name)
                
        

        obj_names_col = Column(processed_names, name='obj_names')
        obj_names_tab = Table([obj_names_col])
        obj_names_tab.write(f'{coadd_directory}processed_obj_names.csv', format='csv', overwrite=True)
                
                
#%% get all gifs 
def get_all_gifs(coadd_directory, obj_names, bands, epochs_tab_path, endswith='_processed.fits', do_sqrt=True, cmap='turbo', scale_values=False):  
           
    for name in obj_names:
        for band in bands:
            print(name)
            
            file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
            # read in obj name and band
            split_path = file_dir.split('/')
            name = split_path[6]
            band = split_path[8]
            endswith_base=endswith.replace('.fits', '')
            if os.path.exists( f'/users/dzakaria/DATA/dzfiles/animations/{name}_{band}_turbo_{endswith_base}.gif'):
                print('gif_exists')
            else:

                make_gif(file_dir, epochs_tab_path=epochs_tab_path, 
                         endswith=endswith, do_sqrt=do_sqrt, cmap=cmap, scale_values=scale_values)
                
#%% make all plots

def make_all_plots(coadd_directory, obj_names, bands):
    
    for name in obj_names:
        for band in bands:
            print(name)
            
            file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
            # read in obj name and band
            split_path = file_dir.split('/')
            name = split_path[6]
            band = split_path[8]
            
            # make the four light cuves
            
            
            """ 
            TARGET AND COMPSTARS
            """
            
            endswith = '_processed.fits'
            endswith_base=endswith.replace('.fits', '')
            if os.path.exists(f'/users/dzakaria/DATA/dzfiles/lightcurves/{name}_{band}{endswith_base}_target_and_compstars_lc.png'):
                print('plot_exists')
            else:
                l, u = make_light_curve(file_dir, endswith=endswith, just_target=False)
            
            
            
            """ 
            JUST THE TARGET (NO COMP STARS) 
            """
            
            endswith = '_processed.fits'
            endswith_base=endswith.replace('.fits', '')
            if os.path.exists( f'/users/dzakaria/DATA/dzfiles/lightcurves/{name}_{band}{endswith_base}_target_lc.png'):
                print('plot_exists')
            else:
                make_light_curve(file_dir, endswith=endswith, just_target=True, lower_value=l, upper_value=u)
                
                
            """ 
            SUBTRACTED IMAGES TARGET LC
            """
            
            endswith = '_subtracted.fits'
            endswith_base=endswith.replace('.fits', '')
            if os.path.exists( f'/users/dzakaria/DATA/dzfiles/lightcurves/{name}_{band}{endswith_base}_target_lc.png'):
                print('plot_exists')
            else:
                make_light_curve(file_dir, endswith=endswith, just_target=True, lower_value=l, upper_value=u)
                
            
            """ 
            DIVIDED IMAGES TARGET LC
            """
            
            endswith = '_subtracted.fits'
            endswith_base=endswith.replace('.fits', '')
            if os.path.exists( f'/users/dzakaria/DATA/dzfiles/lightcurves/{name}_{band}{endswith_base}_target_lc.png'):
                print('plot_exists')
            else:
                make_light_curve(file_dir, endswith=endswith, just_target=True, lower_value=l, upper_value=u)
                
                
#%% fix_ts_tab (with sub and div columns for target phot)

def fix_ts_tab(fits_dir,epochs_tab_path):
    """
    1) open the time series table
    2) display the coadd image with apertures plotted
    3) choose which apertures to keep
    4) choose y min and maxes for plots
    5) save plot
    
    """
  
    # read in obj name and band
    split_path = fits_dir.split('/')
    name = split_path[6]
    band = split_path[8]
    
    time_series_tab_path = f'{fits_dir}/time_series_{name}_{band}.ecsv'
    

    # open time_series_tab
    ts_tab=QTable.read(time_series_tab_path)
    
    colnames = ts_tab.colnames
    
    target_apertures=[]
    compstar_apertures=[]
    compstar_annuli=[]
    for ap_name in colnames:
        if "TARGET" in ap_name:
            ra = ts_tab[ap_name].meta['ra[deg]']
            dec = ts_tab[ap_name].meta['dec[deg]']
            a = ts_tab[ap_name].meta['semi-major_axis[pix]']
            b = ts_tab[ap_name].meta['semi-minor_axis[pix]']
            theta = ts_tab[ap_name].meta['theta[rad]']
            
            aperture_position = SkyCoord(ra, dec, unit='deg')

            aperture = SkyEllipticalAperture(aperture_position, a, b, theta)
            target_apertures.append(aperture)
        elif 'time' not in ap_name: 
            ra = ts_tab[ap_name].meta['ra[deg]']
            dec = ts_tab[ap_name].meta['dec[deg]']
            rad = ts_tab[ap_name].meta['aperture_rad[pix]']
            aperture_position = SkyCoord(ra, dec, unit='deg')
            aperture = SkyCircularAperture(aperture_position, rad)
            compstar_apertures.append(aperture)
            
            
            r_in= ts_tab[ap_name].meta['annulus_inner_rad[pix]']
            r_out= ts_tab[ap_name].meta['annulus_outer_rad[pix]']
            annulus_aperture = SkyCircularAnnulus(aperture_position, r_in=r_in, r_out=r_out)
            compstar_annuli.append(annulus_aperture)
            
    ts_tab_new = time_series_tab(fits_dir, epochs_tab_path, compstar_apertures, compstar_annuli, target_apertures)
    ts_tab_new.write(f'{fits_dir}/time_series_{name}_{band}.ecsv', format='ascii.ecsv', overwrite=True)
    ts_tab_new.write(f'{fits_dir}/time_series_{name}_{band}.csv', format='ascii.csv', overwrite=True)          
    
              