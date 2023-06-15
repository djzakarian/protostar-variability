#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:45:23 2023

# GOAL:
    
    Go through the protostar target list and acquire their single exposure 
    metadata. From this metadata, determine the epochs of each object, 
    and construct the urls for each epoch to download the coadds from IRSA.
    Note: if the coadds on command line timeout, we will need to manually 
    download images using the web interface of the coadder tool (https://irsa.ipac.caltech.edu/applications/ICORE/)
    

@author: dzakaria 
"""

# %% # imports

import pandas as pd
import numpy as np
import pyvo
from astropy.coordinates import SkyCoord 
import astropy.units as u
from astropy.table import vstack, QTable
from astropy.io import ascii
# %%  # Read in the table of targets

# directory path

directory = '/users/dzakaria/dzfiles/'
filename = 'targets-daphne.csv'
path = directory + filename
targets_tab = QTable.read(path, format = 'ascii.csv')

# first column has a weird unicode thing happening - fix that by renaming column
# you can see column names by doing table.colnames()
targets_tab.rename_column('\ufeffobj_name', 'obj_name')

# %%  # Format table (assign units and make SkyCoord object column)

# assign units to columns
targets_tab['ra'].unit= u.hourangle
targets_tab['dec'].unit= u.deg

# make skycoord object column
targets_tab['coord'] = SkyCoord(targets_tab['ra'], targets_tab['dec'], frame='icrs', unit=(u.hourangle, u.deg))



# %% # query_IRSA function 
# make a TAP query to IRSA given: catalog, ra, dec, and radius

def query_IRSA(catalog, ra, dec, rad = 0.5 , columns = '*' ):
    # query base:
    query_base = """
               SELECT {columns}
               FROM {catalog}
               WHERE CONTAINS(POINT('ICRS',crval1, crval2), CIRCLE('ICRS',{ra},{dec},{rad}))=1
               """
    query = query_base.format(columns=columns, catalog=catalog, ra=ra, dec=dec, rad=rad)
    
    # make connection with IRSA
    service = pyvo.dal.TAPService('https://irsa.ipac.caltech.edu/TAP')
    result = service.run_async(query)
    tab = result.to_table()
    return tab


#%% # Make a new output table that will contain the object info as well as epochs

# use the targets_tab to read in the correct columns
epochs_tab = QTable(names=('obj_name', 'ra', 'dec', 'coord', 'obj_epoch', 'date_obs1', 'mjd_obs1', 'date_obs2', 'mjd_obs2', 'n_images'))

epochs_tab['obj_name'].dtype='<U64'
epochs_tab['ra'].dtype='<U64'
epochs_tab['dec'].dtype='<U64'
epochs_tab['date_obs1'].dtype='<U64'
epochs_tab['date_obs2'].dtype='<U64'
#%% # Loop through targets and determine object epochs


# initialize epochs_tab_rownum counter
epochs_tab_rownum = 0

for targets_row in range(len(targets_tab)):
  
    obj_name = targets_tab[targets_row]['obj_name']
    
    
    # print an update every time the next target is processed
    print('target: ', obj_name)
    

    coord = targets_tab[targets_row]['coord']
    
    # query ALLSKY (original WISE observations)
    # WISE All-Sky Single Exposure (L1b) Image Inventory Table (allsky_4band_p1bm_frm)

    allwise_tab = query_IRSA(catalog='allsky_4band_p1bm_frm', ra=coord.ra.deg,
                            dec=coord.dec.deg, rad=0.5, 
                            columns = 'date_obs, mjd_obs, band')
    
    # query NEOWISE (post-cryo WISE observations)
    # query NEOWISE-R Single Exposure (L1b) Image Inventory Table  (wise.wise_allwise_p3am_cdd)

    neowise_tab = query_IRSA(catalog='neowiser_p1bm_frm', ra=coord.ra.deg,
                            dec=coord.dec.deg, rad=0.5,
                            columns = 'date_obs, mjd_obs, band')
 
    epoch = 0
    
    # combine the lists of datses from ALLWISE and NEOWISE and sort by date 
    all_obs_tab = vstack([allwise_tab, neowise_tab])

    all_obs_tab.sort('mjd_obs')
    
    # ===================================================================
    # now: it's time to determine epochs based on the observation dates
    
    # initialize the epoch counter 
    epoch = 0
    
    # add a row to epochs_tab and update rownum counter
    epochs_tab.add_row()
    
    
    # the all_obs_tab is sorted, so naturally, the first observation will be 
    # the first observation of the first epoch
    # before entering the loop, we need to populate the output table with this information
    epochs_tab['date_obs1'][epochs_tab_rownum ] = all_obs_tab[0]['date_obs']
    epochs_tab['mjd_obs1'][epochs_tab_rownum ] = all_obs_tab[0]['mjd_obs']

    epochs_tab['obj_epoch'][epochs_tab_rownum ]=epoch
    epochs_tab['obj_name'][epochs_tab_rownum ] = obj_name
    epochs_tab['ra'][epochs_tab_rownum ] = coord.ra
    epochs_tab['dec'][epochs_tab_rownum ] = coord.dec
    
    # initialize n_images counter 
    # tells us how many images were observed in one epoch (in all bands)
    n_images = 0
    
    for row in range(1, len(all_obs_tab)):
        
        
        
        # compare row to previous row to determine if it is a new epoch or not
        prev_row = row -1
        
        prev_mjd = all_obs_tab[prev_row]['mjd_obs']
        mjd = all_obs_tab[row]['mjd_obs']
        
        
        # when loop reaches the last observation, that should be the final end of the last epoch
        if row == len(all_obs_tab) -1:
            # add to the counter on number of images
            n_images+=1
            
            # first update the end of the previous epoch 
            epochs_tab['date_obs2'][epochs_tab_rownum] = all_obs_tab[prev_row]['date_obs']
            epochs_tab['mjd_obs2'][epochs_tab_rownum] = all_obs_tab[prev_row]['mjd_obs']
            epochs_tab['n_images'][epochs_tab_rownum] = n_images 
            
            
            break
            
            
        
        # if time between obs is less than 2 days... the two obs are in the same epoch
        # we do want to count how many images are in each epoch
        if mjd - prev_mjd <= 15:
            n_images+=1
            continue
        
        # if obs are more than 2 days apart... there's a new epoch!
        elif mjd - prev_mjd > 15:
            
            # add to the counter on number of images
            n_images+=1
            
            # first update the end of the previous epoch 
            epochs_tab['date_obs2'][epochs_tab_rownum] = all_obs_tab[prev_row]['date_obs']
            epochs_tab['mjd_obs2'][epochs_tab_rownum] = all_obs_tab[prev_row]['mjd_obs']
            epochs_tab['n_images'][epochs_tab_rownum] = n_images
            
            # then update the beginning of the next epoch
            epoch +=1 
            
            epochs_tab_rownum += 1
            epochs_tab.add_row()
            
            epochs_tab['date_obs1'][epochs_tab_rownum] = all_obs_tab[row]['date_obs']
            epochs_tab['mjd_obs1'][epochs_tab_rownum] = all_obs_tab[row]['mjd_obs']
            epochs_tab['obj_epoch'][epochs_tab_rownum]=epoch
            epochs_tab['obj_name'][epochs_tab_rownum] = obj_name
            epochs_tab['ra'][epochs_tab_rownum] = coord.ra
            epochs_tab['dec'][epochs_tab_rownum] = coord.dec
            
            
          
            # reset n_images to 0
            n_images = 0
            continue 


#%% save epochs_tab table 
ascii.write(epochs_tab, '{path}15-06-23_epochs_table.ecsv'.format(path=directory), format='ecsv', overwrite=True)
ascii.write(epochs_tab, '{path}15-06-23_epochs_table.csv'.format(path=directory), format='csv', overwrite=True)   


#%% read in tables

epochs_tab=QTable.read('{path}15-06-23_epochs_table.ecsv'.format(path=directory))

#%% clean table

# it had a bunch of empty rows at the end; delete all of those

for row in range(len(epochs_tab)):
    if epochs_tab[row]['n_images']==0:
        epochs_tab.remove_row(row)
        
        
#%% make the urls to coadd for each epoch for bands 1 and 2