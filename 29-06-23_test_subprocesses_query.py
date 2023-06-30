#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 12:43:07 2023

@author: dzakaria
"""

import subprocess

url= "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=20:35:46.450%20+67:53:04.200&band=2&sizeX=0.0667&sizeY=0.0667&date1=15Jul2014%2003:49:41&date2=22Jul2014%2001:36:22&mode=PI"
# url = 'https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=04:39:53.880 +26:03:09.498&band=1&sizeX=0.0667&sizeY=0.0667&date1=22Feb2010 07:21:51&date2=27Feb2010 08:02:59&mode=PI'
file_path = "/users/dzakaria/DATA/dzfiles/TESTING_DOWNLOAD_FUNCTION"
max_attempts = 2

command = ['wget', '--timeout=10', '-O', file_path, url]

try:
    subprocess.run(command, check=True)
    print("File downloaded successfully")


except subprocess.TimeoutExpired:
    print("download took too long. moving on")
    
except subprocess.CalledProcessError as e:
    print(f"failed to download the file. Error: {e}")






#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:20:30 2023

@author: dzakaria
"""

import os

def download_file(url, output_path, max_attempts):
    
    download_sucess=True
    command = f'wget {url} -O {output_path} 2>&1'
    stream = os.popen(command)
    output = stream.read()
    
    if "Saving to:" in output: 
        download_sucess=True
    
    if f"(try:{max_attempts})" in output:
        print(f"reached attempt {max_attempts}. Giving up on download")
        download_sucess=False
    
    return download_sucess
    
    
#%% 
url= "https://irsa.ipac.caltech.edu/cgi-bin/ICORE/nph-icore?locstr=20:35:46.450%20+67:53:04.200&band=2&sizeX=0.0667&sizeY=0.0667&date1=15Jul2014%2003:49:41&date2=22Jul2014%2001:36:22&mode=PI"
output_path = "/users/dzakaria/DATA/dzfiles/TESTING_DOWNLOAD_FUNCTION"
max_attempts = 2

#%%
download_file(url=url, output_path=output_path, max_attempts=max_attempts)
#%%


command = 'wget {url} -O {output_path} 2>&1'.format(url=url, output_path=output_path)
stream = os.popen(command)
output = stream.read()
print(output)

if "Saving to:" in output: 
    download_sucess=True

if "(try:{max_attempts})".format(max_attempts=max_attempts) in output:
    print(f"reached attempt {max_attempts}. Giving up on download")
    download_sucess=False