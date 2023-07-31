
#%% imports
from astropy.table import Table
#%% epochs_tab path


directory ='/users/dzakaria/DATA/dzfiles/'
coadd_directory ='/users/dzakaria/DATA/dzfiles/coadds-0_0667/'



size ='0_0667'

epochs_tab_path='{path}30-06-23_epochs_table_url_size{size}.csv'.format(path=directory, size=size)
#%% list of objects, bands, cmaps

coadd_directory = '/users/dzakaria/DATA/dzfiles/coadds-0_0667/'


obj_names = ['B335','BHR7_IRAS08124-3422','BHR71','CB17','CB230','CB244','CB6',
            'CB68','CG30','Ced110IRS4','CepheusE','DC303.8-14.2','HH111MMS',
            'HH270VLA1','HH46_47','IRAS03282+3035','IRAS03292+3039','IRAS04166+2706',
            'IRAS04169+2702','IRAS04302+2247','IRAS04325+2402','IRAS05295+1247',
            'IRAS05329-0505','IRAS09449-5052','IRAS11072-7727','IRAS15398-3359',
            'IRAS16253-2429','IRAS02086','L1152','L1157','L1165','L1251A','L1251B',
            'L1251C','L1448IRS2','L1448IRS3','L1448-mm','L1489IRS','L1521F',
            'L1527_IRAS04368+2557','L1551IRS5','L1551NE','L1616MMS1A','L1634',
            'L483','L723_IRAS19156+1906','L778','RCrAIRAS32','Serpens1','SerpensMMS3']

# obj_names_tab = Table([obj_names], names=['obj_names'])
# obj_names_tab.write(f'{coadd_directory}obj_names.csv', format='csv', overwrite=True)


# obj_names = ['L1527_IRAS04368+2557']

bands = ['b1', 'b2']

processed_obj_names_tab = Table.read(f'{coadd_directory}obj_names.csv', format='csv')

processed_obj_names_column = processed_obj_names_tab['obj_names']
processed_obj_names  = processed_obj_names_column.tolist()

# come back  'CB68','DC303.8-14.2','HH111MMS','IRAS03292+3039','IRAS04325+2402','IRAS11072-7727', 
#           'L1152','L1251B','L1448-mm', 'L1551IRS5', 'L483','L723_IRAS19156+1906','L778','RCrAIRAS32','Serpens1',
missing = ['BHR7_IRAS08124-3422','L1448IRS3', 'L1448-mm','L778']

# obj_names.remove(processed_obj_names)

# run_pipeline(coadd_directory,missing, bands, epochs_tab_path)
# run_pipeline_compstars_mags(coadd_directory, obj_names, bands, epochs_tab_path)
# get_all_gifs(coadd_directory, missing, bands, epochs_tab_path, endswith='_processed.fits', do_sqrt=False, cmap = 'turbo', scale_values=False)
# make_all_plots(coadd_directory, ['L1527_IRAS04368+2557'], bands)

#%%

for name in obj_names:
    for band in bands:
        print(name)
        file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
        fix_ts_tab(file_dir,epochs_tab_path)
 


#%% only divide pipeline

def run_pipeline_compstars_mags(coadd_directory, obj_names, bands, epochs_tab_path) :
    

    for name in obj_names:
        for band in bands:
            print(name)
            
            file_dir = f'{coadd_directory}{name}/mosaic-int_fits/{band}'
            
            if band == 'b1':
                wcs_apertures, wcs_annuli, target_wcs_apertures = pipeline(file_dir, epochs_tab_path = epochs_tab_path,
                                                     process_files = False, comp_star_phot = True,
                                                    combine = False, subtract = False, divide = False)
                
            else:
                pipeline(file_dir, epochs_tab_path = epochs_tab_path, wcs_apertures = wcs_apertures, 
                         wcs_annuli = wcs_annuli, target_wcs_apertures=target_wcs_apertures,
                         process_files = False, comp_star_phot = True,
                         combine = False, subtract = False, divide = False)
                
                # # after successfully completed, remove the object name from list of files and save to coadd directory
                # obj_names.remove(name)
                # obj_names_col = Column(obj_names, name='obj_names')
                # obj_names_tab = Table([obj_names_col])
                # obj_names_tab.write(f'{coadd_directory}obj_names.csv', format='csv', 
                
                
#%% 
obj_names = ['B335']
                
#%% mag fix pipeline               
run_pipeline_compstars_mags(coadd_directory, missing, bands, epochs_tab_path)

#%%
