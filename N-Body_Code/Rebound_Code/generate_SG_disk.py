# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:32:02 2024

@author: ryanw
"""
import sys
from SGDisk_Model import *

# M = int(input("Input power of SMBH Mass: "))
# f_edd = float(input("Input fraction of Eddington accretion: "))
# alpha = float(input("Input viscosity parameter: "))
Mass = int(sys.argv[1])
Mass = 10.**Mass
f_edd = float(sys.argv[2])
alpha = float(sys.argv[3])
b = 0          # assume that viscosity is proportional to total pressure

disk_params = disk_model(Mass, f_edd, alpha, b)
folder = '/disk_models/'
save_disk_model(disk_params, location=folder, name=f'M{int(np.log10(Mass))}-f{f_edd}-a{alpha}-b{b}')
plot_disk_model(disk_params, save=True, location=folder)

print("Done")