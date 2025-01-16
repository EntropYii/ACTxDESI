import numpy as np
import math
import matplotlib.pyplot as plt
import os
import pickle
#from enlib import enmap


def get_group(i,Ngroups,population):
	lenGroup=len(population)/Ngroups
	sel=np.ones(len(population),dtype=bool)
	if i==0:
		sel[0]=False
	else:
		sel[i*int(lenGroup):(i+1)*int(lenGroup)]=False
	return sel

def estimatorFunction(dt,divsmap,sel):
	return np.sum(np.multiply(dt[sel],divsmap[sel]))/np.sum(divsmap[sel])

with open("Bin_dic.pkl", "rb") as pickle_file:
    bin_dict = pickle.load(pickle_file)


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

#bin_name = config["lum_bin_Jack"][0]
#print(bin_name)

lum_bin = config["lum_bin_Jack"]
dts = []
errs = []
for bin_name in lum_bin:
    if bin_name == "":
        continue
    ras=bin_dict[bin_name][0]
    dec=bin_dict[bin_name][1]
    lum=bin_dict[bin_name][2]
    zs=bin_dict[bin_name][3]
    disk_mean=bin_dict[bin_name][4]
    disk_std=bin_dict[bin_name][5]
    annulus_mean=bin_dict[bin_name][6]
    annulus_std=bin_dict[bin_name][7]
    divs=bin_dict[bin_name][8]

    dt=np.subtract(disk_mean,annulus_mean)


    divsbin=np.array(divs)
    avgs=[estimatorFunction(dt,divsbin,get_group(j,500,dt)) for j in range(len(ras))]
    err=np.sqrt(((len(avgs)-1.)/len(avgs))*np.sum((avgs-(np.sum(np.multiply(dt,divsbin))/np.sum(divsbin)))**2.0))
    dtl=np.sum(np.multiply(dt,divsbin))/np.sum(divsbin)
    
    
    print('-----------------------------')
    print("Luminosity bin:", bin_name)
    print('len(avgs)',len(avgs))
    print('N=',len(zs))
    print('<z>=',np.mean(zs))
    print('dt',dtl )
    print('jk err',err)

    dts.append(dtl)
    errs.append(err)

np.savetxt('AP_dt_err.txt',np.array([dts,errs]))
print('-----------------------------')
print("dt and jk err saved to AP_dt_err.txt")
print('-----------------------------')
