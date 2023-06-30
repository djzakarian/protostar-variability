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

# Show plots in the notebook
# %matplotlib inline

#%% style

# plt.style.use('../photutils_notebook_style.mplstyle')

#%% retrieve data from Hubble

url = 'https://archive.stsci.edu/pub/hlsp/xdf/hlsp_xdf_hst_acswfc-60mas_hudf_f435w_v1_sci.fits'
with fits.open(url) as hdulist:
    hdulist.info()
    data = hdulist[0].data
    header = hdulist[0].header

#%% define the mask
mask = data == 0

#%% store image 
unit = u.electron/u.s
xdf_image=CCDData(data, unit=unit, meta=header, mask=mask)

#%% look at the data
# Set up the figure with subplots
fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))

# Set up the normalization and colormap
norm_image = ImageNormalize(vmin=1e-4, vmax=5e-2, stretch=LogStretch(), clip=False)
cmap = plt.get_cmap('viridis')
cmap.set_over(cmap.colors[-1])
cmap.set_under(cmap.colors[0])
cmap.set_bad('white') # Show masked data as white
xdf_image_clipped = np.clip(xdf_image, 1e-4, None) # clip to plot with logarithmic stretch

# Plot the data
fitsplot = ax1.imshow(np.ma.masked_where(xdf_image.mask, xdf_image_clipped), 
                      norm=norm_image, cmap=cmap)

# Define the colorbar and fix the labels
cbar = plt.colorbar(fitsplot, fraction=0.046, pad=0.04, ticks=LogLocator(subs=range(10)))
labels = ['$10^{-4}$'] + [''] * 8 + ['$10^{-3}$'] + [''] * 8 + ['$10^{-2}$']
cbar.ax.set_yticklabels(labels)

# Define labels
cbar.set_label(r'Flux Count Rate ({})'.format(xdf_image.unit.to_string('latex')), 
               rotation=270, labelpad=30)
ax1.set_xlabel('X (pixels)')
ax1.set_ylabel('Y (pixels)')

#%% creating apertures

""" CREATING APERTURES """

#%% imports

from photutils import find_peaks
from photutils.centroids import centroid_2dg

#%% calculate statistics
mean, median, std = sigma_clipped_stats(xdf_image.data, sigma=3.0, maxiters=5, mask=xdf_image.mask)

#%% find peaks
sources_findpeaks = find_peaks(xdf_image.data, mask=xdf_image.mask, 
                               threshold=20.*std, box_size=30, 
                               centroid_func=centroid_2dg)     
# Display the table
sources_findpeaks

#%% plot the centroids of each of the sources

# Set up the figure with subplots
fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))

# Plot the data
fitsplot = ax1.imshow(np.ma.masked_where(xdf_image.mask, xdf_image_clipped), norm=norm_image)
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

#%% Circular Apertures

""" Circular Apertures """

from photutils import CircularAperture

#%% define the aperture

position = np.column_stack((sources_findpeaks['x_centroid'], sources_findpeaks['y_centroid']))
radius = 10.
circular_aperture = CircularAperture(position, r=radius)


#%%

# Set up the figure with subplots
fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))

# Plot the data
fitsplot = ax1.imshow(np.ma.masked_where(xdf_image.mask, xdf_image_clipped), norm=norm_image)

# Plot the apertures
circular_aperture.plot(color='red', alpha=0.7)

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

# Crop to show an inset of the data
ax1.set_xlim(2000, 3000)
ax1.set_ylim(2000, 1000)


#%% Elliptical Apertures

""" Elliptical Apertures """ 

from photutils import (detect_sources, source_properties, \
                       EllipticalAnnulus, EllipticalAperture)

#%% use centroid, semimajor and semiminor axes and orientation values from source_properties to generate elliptical apertures

# Define threshold and minimum object size
threshold = 5. * std
npixels = 15

# Create a segmentation image
segm = detect_sources(xdf_image.data, threshold, npixels)

# Create a catalog using source properties
catalog = source_properties(xdf_image.data, segm)
table = catalog.to_table()

# Display the table
print(table)

#%% 

r = 3.  # approximate isophotal extent of semimajor axis

# Create the apertures
elliptical_apertures = []
for obj in catalog:
    position = (obj.xcentroid.value, obj.ycentroid.value)
    a = obj.semimajor_axis_sigma.value * r
    b = obj.semiminor_axis_sigma.value * r
    theta = obj.orientation.value
    
    elliptical_apertures.append(EllipticalAperture(position, a, b, theta=theta))
    
#%%

# Set up the figure with subplots
fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))

# Plot the data
fitsplot = ax1.imshow(np.ma.masked_where(xdf_image.mask, xdf_image_clipped), norm=norm_image)

# Plot the apertures
for aperture in elliptical_apertures:
    aperture.plot(color='red', alpha=0.7, axes= ax1)

# Define the colorbar
cbar = plt.colorbar(fitsplot, fraction=0.046, pad=0.04, ticks=LogLocator(subs=range(10)))
labels = ['$10^{-4}$'] + [''] * 8 + ['$10^{-3}$'] + [''] * 8 + ['$10^{-2}$']
cbar.ax.set_yticklabels(labels)

# Define labels
cbar.set_label(r'Flux Count Rate ({})'.format(xdf_image.unit.to_string('latex')), 
               rotation=270, labelpad=30)
ax1.set_xlabel('X (pixels)')
ax1.set_ylabel('Y (pixels)')
ax1.set_title('Elliptical Apertures')

# Crop to show an inset of the data
ax1.set_xlim(2000, 3000)
ax1.set_ylim(2000, 1000)


#%% Sky Coordinates and Apertures

""" Sky Coordinates and Apertures """

from astropy.wcs import WCS

#%% generate sky elliptical apertures

## using the above work included, we could generate sky elliptical apertures in a few lines:
wcs = WCS(header)
sky_elliptical_apertures = [ap.to_sky(wcs) for ap in elliptical_apertures]

from photutils import SkyEllipticalAperture
from astropy.coordinates import SkyCoord

r = 3.  # approximate isophotal extent of semimajor axis

# Create the apertures
sky_elliptical_apertures = []
for obj in catalog:
    # Convert the centroids into RA/Dec using WCS
    ra, dec = wcs.all_pix2world(obj.xcentroid.value, obj.ycentroid.value, 0)
    # Convert the positions to an Astropy SkyCoord object, with units!
    sky_position = SkyCoord(ra, dec, unit=u.deg)
    
    # Define the elliptical parameters, now with units
    a = obj.semimajor_axis_sigma.value * r * u.pix
    b = obj.semiminor_axis_sigma.value * r * u.pix
    theta = obj.orientation.value  * u.rad
    
    # Convert the theta from radians from X axis to the radians from North 
    x_to_north_angle = (90. + header['ORIENTAT']) * u.deg
    x_to_north_angle_rad = x_to_north_angle.to_value(u.rad) * u.rad
    theta -= x_to_north_angle_rad
    
    # Define the apertures
    ap = SkyEllipticalAperture(sky_position, a, b, theta=theta)
    sky_elliptical_apertures.append(ap)

    # note: you can't plot sky apertures

