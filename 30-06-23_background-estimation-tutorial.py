#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:18:49 2023

Background Estimation

@author: dzakaria
"""
#%% make data
from photutils.datasets import make_100gaussians_image

data = make_100gaussians_image() # mean is 5, stdev is 2

#%% plot image

import matplotlib.pyplot as plt

from astropy.visualization import SqrtStretch

from astropy.visualization.mpl_normalize import ImageNormalize

norm = ImageNormalize(stretch=SqrtStretch())

plt.imshow(data, norm=norm, origin='lower', cmap='Greys_r',

           interpolation='nearest')

#%% print median and biweight 
import numpy as np
from astropy.stats import biweight_location
print(np.median(data))  

print(biweight_location(data))
#%% median absolute deviation to est the background noise
from astropy.stats import mad_std
print(mad_std(data))  

#%% Sigma Clipping Sources

""" Sigma Clipping Sources """

#%% bestter est of background
from astropy.stats import sigma_clipped_stats

mean, median, std = sigma_clipped_stats(data, sigma=3.0)

print((mean, median, std))  
#%% masking sources

from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
threshold = detect_threshold(data, nsigma=2.0, sigma_clip=sigma_clip)
segment_img = detect_sources(data, threshold, npixels=10)
footprint = circular_footprint(radius=10)
mask = segment_img.make_source_mask(footprint=footprint)
mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
print((mean, median, std))  

#%% create test image by adding a gradient background to image above
ny, nx = data.shape
y, x = np.mgrid[:ny, :nx]
gradient = x * y / 5000.0
data2 = data + gradient
plt.imshow(data2, norm=norm, origin='lower', cmap='Greys_r',
           interpolation='nearest')  

#%% create a background 2d object using a box size 50x50 and a 3x3 median filter

from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(data2, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

#%% the global median value of the low res background and background rms image 
# can be accessed with the background median and background rms median

print(bkg.background_median)  

print(bkg.background_rms_median)  
#%% plot background image
plt.imshow(bkg.background, origin='lower', cmap='Greys_r',

           interpolation='nearest')

#%% plot background subtracted image

plt.imshow(data2 - bkg.background, norm=norm, origin='lower',

           cmap='Greys_r', interpolation='nearest')
#%% Masking
""" Masking """
#%% create a rotated image that has blank areas and plot it 
from scipy.ndimage import rotate

data3 = rotate(data2, -45.0)

norm = ImageNormalize(stretch=SqrtStretch())  

plt.imshow(data3, origin='lower', cmap='Greys_r', norm=norm,

           interpolation='nearest')  

#%% create acoverage mask and input it into background2d to exclude the regions where we have no data
coverage_mask = (data3 == 0)

bkg3 = Background2D(data3, (15, 15), filter_size=(3, 3),

                    coverage_mask=coverage_mask, fill_value=0.0,

                    exclude_percentile=50.0)

#%% 
norm = ImageNormalize(stretch=SqrtStretch())  

plt.imshow(bkg3.background, origin='lower', cmap='Greys_r', norm=norm,

           interpolation='nearest')

#%% 
norm = ImageNormalize(stretch=SqrtStretch())

plt.imshow(data3 - bkg3.background, origin='lower', cmap='Greys_r',

           norm=norm, interpolation='nearest')


#%% Plotting Meshes
""" Plotting Meshes """
#%%
plt.imshow(data3, origin='lower', cmap='Greys_r', norm=norm,

           interpolation='nearest')

bkg3.plot_meshes(outlines=True, marker='.', color='cyan', alpha=0.3)

plt.xlim(0, 250)

plt.ylim(0, 250)

#%%
























