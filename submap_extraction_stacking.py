import argparse
import numpy as np
import pandas as pd
import h5py
from astropy.io import fits
import pixell 
from pixell import enmap
import os
import matplotlib.pyplot as plt
import pickle

###################### Argument parse and Reading directories from config.txt ##############################
#parser = argparse.ArgumentParser()
#parser.add_argument("-c", "--catalog", dest="catalog", help="cluster catalog to use", metavar="FILE")
#parser.add_argument("-a","--analysis",dest="analysis",help="what analysis to run",metavar="FILE",action='append')
#args = parser.parse_args()
def load_config(config_file="config.txt"):
    config = {}
    with open(config_file, "r") as file:
        for line in file:
            # Skip empty lines and comments
            if line.strip() and not line.startswith("#"):
                # Ensure the line contains an '=' character
                if "=" in line:
                    key, value = line.strip().split("=", 1)
                    # Handle list parsing
                    if "," in value:
                        config[key.strip()] = [item.strip() for item in value.split(",")]
                    else:
                        config[key.strip()] = value.strip()
                else:
                    print(f"Skipping malformed line: {line.strip()}")
    return config

# Load configuration
config = load_config()

catalog_dir = config.get("catalog_dir")
ilc_map_dir = config.get("ilc_map_dir")
submap_dir = config.get("submap_dir")
lum_bin =  config["lum_bin"]
print(lum_bin)

# Ensure output directory exists
#os.makedirs(catalog_dir, exist_ok=True)

print(f"catalog directory is : {catalog_dir}")
print(f"ilc map directory is: {ilc_map_dir}")
print(f"submap_dir directory is: {submap_dir}")



################################### apply flags from header ############################################
catalog = pd.read_csv(catalog_dir, delimiter=',')  # please manually count the rows of flags 
# for full catalog ski row 17 

# Apply the inverse white noise variance cut (45 µK cut for more conservative cuts)
catalog = catalog[(catalog['divcut'] == 2)]

# Apply the galactic plane masking cut (50% for conservativeness)
catalog = catalog[catalog['galcut'] == 2]

#data = data[data['PS'] == 2]
catalog = catalog[catalog['PS15mJy_cut'] == 1]
catalog = catalog[catalog['PS100mJy_cut'] == 1]

catalog = catalog[catalog['S16ILC'] == 1]

################################### submap extraction  ##############################################

# Replace with the actual column names for RA and Dec in your catalog
ra_column = 'ra'
dec_column = 'dec'

# Initialize a dictionary to store coordinates for each bin
coordinates_by_bin = {bin_name: ([], [], [],[], [],[],[],[],[]) for bin_name in ['L43', 'L61', 'L79', 'L98', 'L116', 'L43D', 'L61D', 'L79D', 'L98D']}

# Iterate through each row to handle overlapping bins
for _, row in catalog.iterrows():
    luminosity = row['lum']  # Replace 'Luminosity' with the actual column name
    z = row['z']
    ra = row[ra_column]
    dec = row[dec_column]
    # Check each bin condition and append coordinates to the corresponding bin
    if luminosity > 4.30e10:
        coordinates_by_bin['L43'][0].append(ra)
        coordinates_by_bin['L43'][1].append(dec)
        coordinates_by_bin['L43'][2].append(luminosity)
        coordinates_by_bin['L43'][3].append(z)
    if luminosity > 6.10e10:
        coordinates_by_bin['L61'][0].append(ra)
        coordinates_by_bin['L61'][1].append(dec)
        coordinates_by_bin['L61'][2].append(luminosity)
        coordinates_by_bin['L61'][3].append(z)
    if luminosity > 7.90e10:
        coordinates_by_bin['L79'][0].append(ra)
        coordinates_by_bin['L79'][1].append(dec)
        coordinates_by_bin['L79'][2].append(luminosity)
        coordinates_by_bin['L79'][3].append(z)
    if luminosity > 9.80e10:
        coordinates_by_bin['L98'][0].append(ra)
        coordinates_by_bin['L98'][1].append(dec)
        coordinates_by_bin['L98'][2].append(luminosity)
        coordinates_by_bin['L98'][3].append(z)
    if luminosity > 11.60e10:
        coordinates_by_bin['L116'][0].append(ra)
        coordinates_by_bin['L116'][1].append(dec)
        coordinates_by_bin['L116'][2].append(luminosity)
        coordinates_by_bin['L116'][3].append(z)
    if 4.30e10 < luminosity <= 6.10e10:
        coordinates_by_bin['L43D'][0].append(ra)
        coordinates_by_bin['L43D'][1].append(dec)
        coordinates_by_bin['L43D'][2].append(luminosity)
        coordinates_by_bin['L43D'][3].append(z)
    if 6.10e10 < luminosity <= 7.90e10:
        coordinates_by_bin['L61D'][0].append(ra)
        coordinates_by_bin['L61D'][1].append(dec)
        coordinates_by_bin['L61D'][2].append(luminosity)
        coordinates_by_bin['L61D'][3].append(z)
    if 7.90e10 < luminosity <= 9.80e10:
        coordinates_by_bin['L79D'][0].append(ra)
        coordinates_by_bin['L79D'][1].append(dec)
        coordinates_by_bin['L79D'][2].append(luminosity)
        coordinates_by_bin['L79D'][3].append(z)
    if 9.80e10 < luminosity <= 11.60e10:
        coordinates_by_bin['L98D'][0].append(ra)
        coordinates_by_bin['L98D'][1].append(dec)
        coordinates_by_bin['L98D'][2].append(luminosity)
        coordinates_by_bin['L98D'][3].append(z)



# Convert lists to arrays for each bin
#for bin_name in coordinates_by_bin:
 #   coordinates_by_bin[bin_name] = (np.array(coordinates_by_bin[bin_name][0]), np.array(coordinates_by_bin[bin_name][1]),np.array(coordinates_by_bin[bin_name][2]),np.array(coordinates_by_bin[bin_name][3])
  #                                  ,np.array(coordinates_by_bin[bin_name][4]), np.array(coordinates_by_bin[bin_name][5]),np.array(coordinates_by_bin[bin_name][6]),np.array(coordinates_by_bin[bin_name][7])
   #                                 ,np.array(coordinates_by_bin[bin_name][8]))



# Check if bin is healthy 
#example_bin = "L43D"  # Replace with the bin you want
#if example_bin in coordinates_by_bin:
    #ra_coords, dec_coords, luminosity, redshift = coordinates_by_bin[example_bin]
    #print(f"Coordinates for bin {example_bin}:")
    #print("RA:", ra_coords)
    #print("Dec:", dec_coords)
    #print("Luminosity:", luminosity)
    #print("redshift:", redshift)
#else:
    #print(f"No data for bin {example_bin}")
Dsum = len(coordinates_by_bin['L43D'][1])+ len(coordinates_by_bin['L61D'][1])+ len(coordinates_by_bin['L79D'][1])+ len(coordinates_by_bin['L98D'][1])
Jdiff = len(coordinates_by_bin['L43'][1])-len(coordinates_by_bin['L116'][1])
if Dsum-Jdiff ==0:
    print('~Bin size looks good ! ')
##################################### define pixel size #############################################
pix_arcmin = 0.5  #each pixel is 0.5 arcmins 
resolution_factor =  1/pix_arcmin # for making sure we get the exact size of submaps we want


##################################### functions for AP masks ########################################
def circular_mask(h, w,center=None, radius=None):
    if center is None:  # use the middle of the image
        center = [int(h/2), int(w/2)]
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], h-center[0], w-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[1])**2 + (Y - center[0])**2)

    mask = dist_from_center <= radius
    return mask

def ring_mask(h, w,center=None, radius_inner=None, radius_outter = None ):
    if center is None:  # use the middle of the image
        center = [int(h/2), int(w/2)]

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[1])**2 + (Y - center[0])**2)

    mask = (dist_from_center <= radius_outter) & (dist_from_center >= radius_inner)
    return mask
##################################### AP params #####################################################

fout = np.sqrt(2) #factor which defines outer AP radius
fres = 2  # resolution factor for inverse of pixel armin representation 
fres = float(fres)

radius = 2.1*fres#arcmin 


radius = float(radius)
radius_out = radius*fout #defining AP radii 

################################### submap extraction  ##############################################

def extract_submaps(fits_file, ra, dec, output_dir, lum_bin, divmap, submap_size=18.0,):

    map_data = enmap.read_map(fits_file, hdu=0)
    wcs = map_data.wcs
    if len(map_data) <= 5: 
        map_data = map_data[0]
    print(len(map_data))
    #catalog_hdulist = fits.open(catalog_file)
    #catalog_data = catalog_hdulist[1].data

    #ra = catalog_data['RADeg']
    #dec = catalog_data['decDeg']
    #mass = catalog_data['M500c']
    
    # below are coordinated copied from table 7 of Hassefield paper
    #ra = [2.0418, 3.0152, 3.7276, 4.4138, 4.5623, 5.5553, 6.5699, 11.1076, 11.3051, 12.7875, 14.5189, 14.7855, 16.2195, 19.9971, 21.8227, 24.8407, 28.1764, 29.1008, 31.5567, 33.8699, 34.5626, 34.9533, 34.9759, 35.3925, 35.7939, 37.1250, 37.7273, 39.9718, 40.0102, 40.3129, 41.4645, 42.5370, 44.1354, 45.2925, 45.4158, 45.8343, 47.0481, 50.1239, 51.7075, 54.2438, 55.5008, 55.6845, 57.1612, 57.1605, 306.3006, 312.6264, 312.6814, 312.7935, 312.7885, 313.8581, 314.7234, 322.1036, 322.4186, 322.5367, 323.7907, 323.8151, 323.9310, 328.2375, 328.6319, 329.0407, 335.1922, 337.3042, 343.3432, 345.6427, 346.9176, 351.8660, 354.4156, 357.9349]
    #dec = [2.0204, -0.7693, -0.9502, -0.8580, -0.3795, -0.6050, 1.3367, 1.2221, -1.8827, 0.9323, 0.5106, -0.8326, 0.0495, 0.9193, 0.3468, -1.4769, 1.0059, -1.3879, -1.2428, 0.5091, -0.6883, 0.3755, 1.4973, -0.2063, -0.9466, 0.5033, -0.4043, -1.5758, 1.2693, -0.3109, -0.7013, 0.1403, 0.1049, -1.1716, 1.9219, 1.9214, 1.0607, 0.5399, -0.7312, -1.1705, 1.0873, -0.2899, 0.4892, -0.4681, 0.5130, -0.9311, 1.3857, 0.9488, 2.2628, 1.0985, 1.3836, 1.5996, 0.0891, 0.7590, -1.0396, 1.4247, 0.1568, -1.2458, -0.8197, 1.3857, -0.7095, -0.0743, -0.5280, 0.0419, 1.5161, -2.0777, 0.2690, 0.1538]
    
    # Define submap size in arcminutes and convert to radians
    submap_size_rad = np.deg2rad(submap_size / 60.0)
    total_items = len(ra)
    milestones = [total_items * i // 10 for i in range(1, 11)]
    for i, (ra_source, dec_source) in enumerate(zip(ra, dec)):

        pos = np.deg2rad([dec_source, ra_source])  # [DEC, RA] in radians
        # bounding box for the submap
        box = np.array([[pos[0] - submap_size_rad / 2, pos[1] - submap_size_rad / 2],
                        [pos[0] + submap_size_rad / 2, pos[1] + submap_size_rad / 2]])
        
        # Check if the current iteration hits a milestone
        if i + 1 in milestones:  # i + 1 because progress is 1-based
            progress = (i + 1) / total_items * 100
            print(f"extraction: {int(progress)}% complete")
            
        submap = enmap.submap(map_data, box=box)
        if submap.shape[0] != 36 or submap.shape[1] != 36: 
            print([dec_source, ra_source], submap.shape)

        h, w = submap.shape
        if divmap == True: 
            mask = circular_mask(h, w, radius=radius)
            div_mean = np.mean(submap[mask])
            coordinates_by_bin[lum_bin][8].append(div_mean)
        else: 
            mask = circular_mask(h, w, radius=radius)
            mask_outter_ring = ring_mask(h, w, radius_inner=radius, radius_outter = radius_out)
            disk_mean = np.mean(submap[mask])
            disk_std = np.std(submap[mask])
            coordinates_by_bin[lum_bin][4].append(disk_mean)
            coordinates_by_bin[lum_bin][5].append(disk_std)

            
            ring_mean = np.mean(submap[mask_outter_ring])
            ring_std = np.std(submap[mask_outter_ring])
            coordinates_by_bin[lum_bin][6].append(ring_mean)
            coordinates_by_bin[lum_bin][7].append(ring_std)

            #tsz_signal = tsz_signal_inner - tsz_signal_outter_ring
            
        
        # Save the submap to a FITS file
        submap_filename = f"{output_dir}/submap_{i}.fits"
        enmap.write_map(submap_filename, submap)
    
    print("~Submap Extraction complete!")

def stack_submaps(submap_dir, output_file):
    # List all submap 
    submap_files = [os.path.join(submap_dir, f) for f in os.listdir(submap_dir) if f.endswith('.fits')]
    
    #get shape and initialize the stack
    first_submap = enmap.read_map(submap_files[0])
    stack = np.zeros_like(first_submap)
    count = 0
    
    bad_count = 0  # for weried shaped submap count 
    total_items = len(submap_files)
    milestones = [total_items * i // 10 for i in range(1, 11)]
    for submap_file in submap_files:
        submap = enmap.read_map(submap_file)
        if submap.shape[0] != 36 or submap.shape[1] != 36 : 
            bad_count += 1
            print("This file", submap_file,"has shape", submap.shape,"Total bad submap", bad_count)
            continue
        if count + 1 in milestones:  # i + 1 because progress is 1-based
            progress = (count + 1) / total_items * 100
            print(f"stacking: {int(progress)}% complete")
        stack += submap
        count += 1
    print("Stacking complete!")
    # Average 
    stack /= count
    
    enmap.write_map(output_file, stack)
    print(f"Stacking completed. Stacked map saved as {output_file}")




resolution_factor = 2
# for Bin_name in ['L116','L98', 'L79','L61', 'L43', 'L43D', 'L61D', 'L79D', 'L98D']:
for Bin_name in lum_bin:
    if Bin_name == "":
        continue
    #fits_file = filename = Home + '/ACTxDESI/s22_product/act_planck_s08_s22_ftot_night_map.fits'
    #catalog_file = '/Users/yi/Documents/CMB_SZ/DR5_cluster-catalog_v1.1.fits'
    
    #base_output_path = "/Users/yi/Documents/CMB_SZ/ACTxDESI/s22_submap/"

    # Combine them to create the full path
    output_dir = f"{submap_dir}{Bin_name}"

    #output_dir = "/Users/yi/Documents/CMB_SZ/ACTxDESI/DR15_submap/L98D/"

    ###################### submaps for CMB #################
    ra = coordinates_by_bin[Bin_name][0]
    dec = coordinates_by_bin[Bin_name][1]
    Lum = coordinates_by_bin[Bin_name][2]
    z = coordinates_by_bin[Bin_name][3]
    print('-----------------------------')
    print('Now start to extract submaps for Lum bin', Bin_name, ':')
    extract_submaps(ilc_map_dir, ra, dec, output_dir, Bin_name, divmap= False)
    Lum_mean = sum(Lum)/len(Lum)
    z_mean = sum(z)/len(z)
    print ('the avergae luminosity and redshift of ', Bin_name,' is',Lum_mean, z_mean,'Total number of sources in this bin is', len(Lum))
        
    ################ submap for inverse variance map #####################
    print('-----------------------------')
    print('Now start to extract Div-submaps for Lum bin', Bin_name, ':')
    extract_submaps(config.get('ilc_div_map_dir'), ra, dec, config.get('div_submap_dir'), Bin_name, divmap= True)

    ######################### naive stacking ###################
    base_path_stacking = config.get("stacked_map_dir")
    output_file =f"{base_path_stacking}{Bin_name}_stacked.fits"
    stack_submaps(output_dir, output_file)  # output_dir here is where submap stored
    
    stacked_map = enmap.read_map(output_file)
    plt.imshow(stacked_map,cmap='viridis')
    ax = plt.gca()
    ax.set_xticks(np.arange(0, stacked_map.shape[1], resolution_factor))
    ax.set_yticks(np.arange(0, stacked_map.shape[0], resolution_factor))

    # Adjust tick labels to count every 2 pixels as one unit
    ax.set_xticklabels(np.arange(0, stacked_map.shape[1] // 2))
    ax.set_yticklabels(np.arange(0, stacked_map.shape[0] // 2))
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('Intensity')
    plt.title(Bin_name)
    plt.xlabel('Arcminutes')
    plt.ylabel('Arcminutes')
    plt.savefig(f"{base_path_stacking}{Bin_name}_stacked.png")

# Export dictionary to a pickle file
print('-----------------------------')
with open("Bin_dic.pkl", "wb") as pickle_file:
    pickle.dump(coordinates_by_bin, pickle_file)
print("Dictionary exported as Pickle!!")
print('-----------------------------')
