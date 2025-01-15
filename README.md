# ACTxDESI Pipeline

This repository contains tools for processing ACT_ILC maps and performing submap extraction, stacking, and jackknife uncertainty estimation. Follow the steps below to set up and run the pipeline.

---

## 1. Download the ACT_ILC Map
Please download the ACT_ILC map from [here](URL). Set up the directory path in the relevant sections of the code.

---

## 2. Catalog and Point Source Mask
Ensure you have the catalog and point source mask produced by Yulin. These will be required for the pipeline.

---

## 3. Submap Extraction and Stacking
To perform submap extraction and stacking:

1. Place all file directories in the `config.txt` file. This file allows you to:
   - Toggle luminosity bins.
   - Specify the catalog file.
   Detailed instructions can be found within the `config.txt` file.

2. Once `config.txt` is set up, run the `submap_extraction_stacking.py` script. This will:
   - Generate submaps of the CMB with accompanying divmaps.
   - Output a stacked image.
   - Create a dictionary containing the following 8 lists for each luminosity bin:
     - `ra`, `dec`, `lum`, `z`
     - `disk_mean`, `disk_std`
     - `ring_mean`, `ring_std`
     - `divsubmap_mean`

### Command to Run:
```bash
python submap_extraction_stacking.py
