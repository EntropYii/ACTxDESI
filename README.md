ACTxDESI pipeline

Please first dowload ACT_ILC map from here:
And set up the directory in the following code.
1) Catalog and point source mask produced by Yulin.


2) Before extracting the submap, please first put all the file directories in the config.txt where you can also toggle the luminosity bins and catalog file, please see the detailed instrution in the file. After setup the config.txt, you can run the submap_extraction_stacking.py which will generate the submap of CMB and accompany divmap and then output a stacking image and also a dictionary that contains 8 lists for each luminosity bins which are ra, dec, lum , z,disk_mean, disk_std,ring_mean, ring_std, divsubmap_mean respectively.  

Run the code, simply in terminal type: python submap_extraction_stacking.py

3) Also before

