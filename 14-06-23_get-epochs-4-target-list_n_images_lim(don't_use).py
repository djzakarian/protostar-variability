#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:45:23 2023

14-06-23_get-epochs-4-target-list.py

# GOAL:
    
    Go through the protostar target list and acquire their single exposure 
    metadata. From this metadata, determine the epochs of each object, 
    and construct the urls for each epoch to download the coadds from IRSA.
    Note: if the coadds on command line timeout, we will need to manually 
    download images using the web interface of the coadder tool (https://irsa.ipac.caltech.edu/applications/ICORE/)
    

    SPLIT epochs if n_images>100
    THIS IS NOT THE FINAL CODE, WE DON"T NEED TO SPLIT EPOCHS HERE
    
@author: dzakaria 
"""

# %% # imports

import numpy as np
import pyvo
from astropy.coordinates import SkyCoord 
import astropy.units as u
from astropy.table import vstack, QTable
from astropy.io import ascii

from datetime import datetime
# %%  # Read in the table of targets

# directory path

directory = '/users/dzakaria/DATA/dzfiles/'
filename = 'targets-daphne.csv'
path = directory + filename
targets_tab = QTable.read(path, format = 'ascii.csv')

# first column has a weird unicode thing happening - fix that by renaming column
# you can see column names by doing table.colnames()
targets_tab.rename_column('\ufeffobj_name', 'obj_name')



#%% some of the targets have imprecise coordinates... fix those here

# L483
targets_tab[2]['ra'] = '18 17 29.8'
targets_tab[2]['dec'] = '-4 39 36.4'

# %%  # Format table (assign units and make SkyCoord object column - also fix names)

# assign units to columns
targets_tab['ra'].unit= u.hourangle
targets_tab['dec'].unit= u.deg

# make skycoord object column
targets_tab['coord'] = SkyCoord(targets_tab['ra'], targets_tab['dec'], frame='icrs', unit=(u.hourangle, u.deg))

# remove slashes and spaces
for i in range(len(targets_tab)):
    targets_tab['obj_name'][i] = targets_tab['obj_name'][i].replace('/','_')
    targets_tab['obj_name'][i] = targets_tab['obj_name'][i].replace(' ', '')


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

        if targets_row >= 7 and n_images >= 100:
            
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
        
        elif mjd - prev_mjd <= 15:
            n_images+=1
            continue
        
        
        # for rows greater than 7 (where my code broke), divide rows with more than 100 images
        # if obs are more than 2 days apart ... there's a new epoch!
        
        
        
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
ascii.write(epochs_tab, '{path}29-06-23_epochs_table.ecsv'.format(path=directory), format='ecsv', overwrite=True)
ascii.write(epochs_tab, '{path}29-06-23_epochs_table.csv'.format(path=directory), format='csv', overwrite=True)   


#%% read in tables

epochs_tab=QTable.read('{path}29-06-23_epochs_table.ecsv'.format(path=directory))

#%% clean table

# it had a bunch of empty rows at the end; delete all of those

for row in reversed(range(len(epochs_tab))):
    if epochs_tab[row]['n_images']==0:
        epochs_tab.remove_row(row)
        
        
#%% in order to save the urls in the table: add new columns for urls in band1 and band2

url_band1 = np.empty((1,), dtype='<U256')
url_band2 = np.empty((1,), dtype='<U256')

epochs_tab['url_band1'] = url_band1
epochs_tab['url_band2'] = url_band2

#%% make the urls to coadd for each epoch for bands 1 and 2

size = 0.0667

for row in range(len(epochs_tab)):
    
    if row%10==0:
        print(row)
    
    ra = epochs_tab[row]['ra']
    dec = epochs_tab[row]['dec']
    date1_str = epochs_tab[row]['date_obs1']
    date2_str = epochs_tab[row]['date_obs2']
    
    # formate date1 and date2 for the url
    
    # date1:
    dt1 = datetime.strptime(date1_str, '%Y-%m-%d %H:%M:%S.%f')
    date1 = dt1.strftime('%d%b%Y %H:%M:%S')
    
    dt2 = datetime.strptime(date2_str, '%Y-%m-%d %H:%M:%S.%f')
    date2 = dt2.strftime('%d%b%Y %H:%M:%S')
    
    coordinate = SkyCoord(ra=ra, dec=dec, unit=u.deg)
    
    coord = coordinate.to_string('hmsdms', sep=':', precision = 3)
        
    
    url_band1 = "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr={coord}&band={band}&sizeX={sizeX}&sizeY={sizeY}&date1={date1}&date2={date2}&mode=PI" \
            .format(coord=coord,band=1, sizeX=size, sizeY=size, date1=date1, date2=date2)
            
    url_band2 = "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr={coord}&band={band}&sizeX={sizeX}&sizeY={sizeY}&date1={date1}&date2={date2}&mode=PI" \
            .format(coord=coord,band=2, sizeX=size, sizeY=size, date1=date1, date2=date2)
            
    
    # update the table with the urls for bands 1 and 2
    epochs_tab['url_band1'][row] = url_band1
    epochs_tab['url_band2'][row] = url_band2
    
    


#%% delete the extra rows

for row in reversed(range(len(epochs_tab))):
    if epochs_tab[row]['mjd_obs1'] == 0:
        epochs_tab.remove_row(row)
    
#%% save files
    
    ascii.write(epochs_tab, '{path}29-06-23_epochs_table_url_size0_0667.ecsv'.format(path=directory), format='ecsv', overwrite=True)
    ascii.write(epochs_tab, '{path}29-06-23_epochs_table_url_size0_0667.csv'.format(path=directory), format='csv', overwrite=True)   

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    