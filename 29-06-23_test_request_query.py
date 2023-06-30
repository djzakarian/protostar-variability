#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:52:14 2023

@author: dzakaria
"""

import requests


url= "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=20:35:46.450%20+67:53:04.200&band=2&sizeX=0.0667&sizeY=0.0667&date1=15Jul2014%2003:49:41&date2=22Jul2014%2001:36:22&mode=PI"
# url = 'https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=04:39:53.880 +26:03:09.498&band=1&sizeX=0.0667&sizeY=0.0667&date1=22Feb2010 07:21:51&date2=27Feb2010 08:02:59&mode=PI'
file_path = "/users/dzakaria/DATA/dzfiles/TESTING_DOWNLOAD_FUNCTION"
max_attempts = 2

# response = requests.get(url, timeout=60)

try:
    
    response = requests.get(url, timeout=600, allow_redirects=True)
    
    if response.status_code == 200:
        with open(file_path, "wb") as file:
            file.write(response.content)
            print('file downloaded successfully')
    
    else: 
        print('failed to download the file')
        
except requests.exceptions.Timeout:
    print("Download took longer than 10 minutes. Process terminated")
    
except requests.exceptions.RequestException as e:

    print('another error occured: ', e )