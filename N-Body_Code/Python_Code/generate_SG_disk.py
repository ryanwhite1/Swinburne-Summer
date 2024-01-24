# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:32:02 2024

@author: ryanw
"""

from SGDisk_Model import *

M = int(input("Input power of SMBH Mass: "))
f_edd = float(input("Input fraction of Eddington accretion: "))
alpha = float(input("Input viscosity parameter: "))
b = 0.          # assume that viscosity is proportional to total pressure

disk_params = disk_model(M, f_edd, alpha, b)
folder = '/disk_models/'
save_disk_model(disk_params, location=folder, name=f'M{int(np.log10(M))}-f{f_edd}-a{alpha}-b{b}')

print("Done")