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








#%%
















