import argparse
import numpy as np
import pandas as pd
import h5py
from astropy.io import fits
import pixell 
from pixell import enmap
import os
import matplotlib.pyplot as plt

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


###################### Import Lum Bin dictionary #############################

with open("coordinates_by_bin.pkl", "rb") as pickle_file:
    loaded_dict = pickle.load(pickle_file)

