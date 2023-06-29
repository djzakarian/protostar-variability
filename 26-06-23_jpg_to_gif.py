#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:49:30 2023

GOAL: 
    For a given object, make GIFs of the images in both bands

@author: dzakaria
"""


#%% imports
import glob
from PIL import Image

#%% directory of images
directory_b1 = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/L1527_IRAS04368+2557/mosaic-int_jpg/b1/*.jpg'
directory_b2 = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/L1527_IRAS04368+2557/mosaic-int_jpg/b1/*.jpg'

#%% create a list of image filenames
image_files_b1 = glob.glob(directory_b1)
image_files_b2 = glob.glob(directory_b2)

#%% create a list to store the frames of the GIF

frames_b1 = []
frames_b2 = []

#%% open each file and append it to the frames list

for file in image_files_b1:
    frame = Image.open(file)
    frames_b1.append(frame)
    
   
for file in image_files_b2:
    frame = Image.open(file)
    frames_b2.append(frame) 
    
# save the frames as an animated GIF
frames_b1[0].save("b1_animated.gif", format='GIF', append_images=frames_b1[1:], save_all=True, duration=400, loop=0)
