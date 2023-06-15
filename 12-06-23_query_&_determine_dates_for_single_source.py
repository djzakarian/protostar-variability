#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 15:17:55 2023


GOAL: 

    Acquire single exposure metadata from WISE and NEOWISE observations. 
    This can be used to determine appropriate windows for coadding images 
    from a given observation epoch. This code finds the appropriate windows
    for a single coordinate. Future code will utilize these methods to loop 
    through a list of coordinates.

@author: dzakaria
"""


#%%

import pyvo # TAP queries from IRSA

import numpy as np


from astropy.coordinates import SkyCoord 
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.table import vstack, Table, QTable


#%%


# define coordinate using SkyCoord

L1527_coord = SkyCoord(ra='04 39 53.97816' , dec = '+26 03 09.6804', unit=(u.hour, u.deg))



#%%

# define function to query IRSA using pyvo
# adjustable parameters: 
#    columns, catalog, coordinate of interest (decimal deg), radius around point (deg)


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


#%%

# query: WISE All-Sky Single Exposure (L1b) Image Inventory Table (allsky_4band_p1bm_frm)

allwise_tab = query_IRSA(catalog='allsky_4band_p1bm_frm', ra=L1527_coord.ra.deg,
                        dec=L1527_coord.dec.deg, rad=0.5,
                        columns = 'date_obs, mjd_obs, band')


# query NEOWISE-R Single Exposure (L1b) Image Inventory Table  (wise.wise_allwise_p3am_cdd)

neowise_tab = query_IRSA(catalog='neowiser_p1bm_frm', ra=L1527_coord.ra.deg,
                        dec=L1527_coord.dec.deg, rad=0.5,
                        columns = 'date_obs, mjd_obs, band')




#%%


# combine the lists of datses from ALLWISE and NEOWISE and sort by date

all_obs_tab = vstack([allwise_tab, neowise_tab])

all_obs_tab.sort('mjd_obs')

#%%

# new astropy table to record the date ranges that are needed for the object
coadd_epochs_tab = QTable(names=('obj_name', 'crval1', 'crval2', 'obj_epoch', 
                             'date_obs1', 'mjd_obs1', 'date_obs2', 'mjd_obs2', 'n_images'))


# change data type of 'obj_name', and 'date_obs' 1&2 columns to a string data type
coadd_epochs_tab['obj_name'].dtype = '<U64'
coadd_epochs_tab['date_obs1'].dtype = '<U64'
coadd_epochs_tab['date_obs2'].dtype = '<U64'

#%%

# initialize the 'epoch' counter to track which epoch we're on in the loop
epoch = 0


# object specific info... will vary when looping through the list of targets
obj_name  = 'L1527'
ra = L1527_coord.ra.deg
dec = L1527_coord.dec.deg


coadd_epochs_tab.add_row()

# the table is sorted, so the first row is the earliest observation
# this is the beginning of the first epoch
coadd_epochs_tab['date_obs1'][0] = all_obs_tab[0]['date_obs']
coadd_epochs_tab['mjd_obs1'][0] = all_obs_tab[0]['mjd_obs']

coadd_epochs_tab['obj_epoch'][0]=epoch
coadd_epochs_tab['obj_name'][0] = obj_name
coadd_epochs_tab['crval1'][0] = ra
coadd_epochs_tab['crval2'][0] =dec
                         
# initilize number of images counter per epoch (this includes all bands!)
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
        coadd_epochs_tab['date_obs2'][-1] = all_obs_tab[prev_row]['date_obs']
        coadd_epochs_tab['mjd_obs2'][-1] = all_obs_tab[prev_row]['mjd_obs']
        coadd_epochs_tab['n_images'][-1] = n_images
        
        
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
        coadd_epochs_tab['date_obs2'][-1] = all_obs_tab[prev_row]['date_obs']
        coadd_epochs_tab['mjd_obs2'][-1] = all_obs_tab[prev_row]['mjd_obs']
        coadd_epochs_tab['n_images'][-1] = n_images
        
        # then update the beginning of the next epoch
        epoch +=1 
        
        coadd_epochs_tab.add_row()
        
        coadd_epochs_tab['date_obs1'][-1] = all_obs_tab[row]['date_obs']
        coadd_epochs_tab['mjd_obs1'][-1] = all_obs_tab[row]['mjd_obs']
        coadd_epochs_tab['obj_epoch'][-1]=epoch
        coadd_epochs_tab['obj_name'][-1] = obj_name
        coadd_epochs_tab['crval1'][-1] = ra
        coadd_epochs_tab['crval2'][-1] = dec
        
        
      
        # reset n_images to 0
        n_images = 0
        continue 
    

#%%

# generate the url that will eventually go to command line

coord = '043953.59 +260305.50'
size = 0.1
date1 = '22Feb2010+01:00:00'
date2 = '27Feb2010+11:00:00'
band = 1

url = "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr={coord} \
    &band={band}&sizeX={sizeX}&sizeY={sizeY}&date1={date1}&date{date2}&mode=PI" \
        .format(coord=coord,band=band, sizeX=size, sizeY=size, date1=date1, date2=date2)


#%%

# add another column for the duration of the epoch

coadd_epochs_tab['duration'] = coadd_epochs_tab['mjd_obs2'] - coadd_epochs_tab['mjd_obs1']

## curl -o out.xml "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=L1527&band=1&date1=22Feb2010+01:00:00&date2=27Feb2010+11:00:00&mode=PI"

## wget -O L1527.xml "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=043953.59 +260305.50&band=1&sizeX=0.05&sizeY=0.05&date1=22Feb2010+01:00:00&date2=27Feb2010+11:00:00&mode=PI"

## wget -O L1527.xml "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=043953.59 +260305.50&band=1&sizeX=0.07&sizeY=0.07&date1=22Feb2010+01:00:00&date2=27Feb2010+11:00:00&mode=PI"



