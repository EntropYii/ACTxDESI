# ACTxDESI Pipeline

## Getting Started

### Step 1: Prepare the Catalog and Point Source Mask
The catalog and point source mask are produced by Yulin. Ensure these are available before proceeding.

### Step 2: Configure and Extract Submaps
1. Before extracting the submap, set up the required file directories in `config.txt`. This file also allows you to:
   - Toggle luminosity bins.
   - Specify the catalog file.
   
   Detailed instructions can be found within `config.txt`.

2. Once `config.txt` is configured, run the submap extraction script by executing:

   ```bash
   python submap_extraction_stacking.py

This script will:

Generate submaps of the CMB along with accompanying divmap.
Output:
A stacking image.
A dictionary containing 8 lists for each luminosity bin: ra, dec, lum, z, disk_mean, disk_std, ring_mean, ring_std, and divsubmap_mean.

### Step 3: Jackknife uncertainty estimates

Update the following parameters in config.txt:
Arcminutes per pixel.
Luminosity bin for uncertainty estimates.
You can enable all_bins = True to run the estimates for all luminosity bins.
Run the jackknife uncertainty estimation script by executing:

    ```bash
    python Jackknife_uncertainty.py

