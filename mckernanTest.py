# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 10:22:15 2023

@author: ryanw
"""

import numpy as np
import matplotlib.pyplot as plt

class AGNDisk(object):
    def __init__(self, smbhmass=1e8):
        self.mass = smbhmass
        self.mass_ratio = smbhmass / 1e8
        self.r_g = 4.3 * 10**-3 * self.mass / (9 * 10**10) # GM/c^2 in units of pc
    
    
    def migration_ts(self, mass, radius):
        ''''''
        surf_dens_ratio = 10**5 / self.disk_surfdens(radius) # (sigma / 10^5)^-1
        mass_ratio = 5 / mass # (mass / 5M_odot)^-1
        scale_height_ratio = (self.disk_aspectratio(radius) / 0.02)**2 # assume center of mass radius is just orbital radius, i.e. M_smbh >> M_bh
        orbit_ratio = (radius / (10**4 * self.r_g))**(-1/2)
        ts = 38 * orbit_ratio * mass_ratio * scale_height_ratio * surf_dens_ratio * self.mass_ratio**(3/2)
        return ts
    
    
    def disk_temp(self, radius):
        ''''''
        logr = np.log10(radius)
        if logr > 0:
            return 10**3.8
        else:
            return 10**(-0.6 * logr + 3.8)
        
        
    def disk_surfdens(self, radius):
        '''+1 in the power to go from g/cm^2 to kg/m^2'''
        logr = np.log10(radius)
        if logr <= -2:
            return 10**(logr + 8 + 1)
        else:
            return 10**(-2 * logr + 2 + 1)
        
        
    def disk_aspectratio(self, radius):
        '''Aspect ratio of scale height: h = H / r'''
        logr = np.log10(radius)
        if logr <= -2:
            return 10**(-0.7 * logr - 3.6)
        else:
            return 10**(0.6 * logr - 1)
        
        
    def disk_optdepth(self, radius):
        ''''''
        logr = np.log10(radius)
        if logr <= -2:
            return 10**(logr + 7)
        else:
            return 10**(-3.5 * logr - 2)
        
        
    def disk_toomre(self, radius):
        ''''''
        logr = np.log10(radius)
        if logr >= -2:
            return 1 
        else:
            return 10**(-5 * logr - 10)

agn = AGNDisk()
def prob_encounter(factor, m1, m2, N_m1, N_m2, r1, r2):
    return factor * N_m1 * N_m2 / (agn.migration_ts(m1, r1) * agn.migration_ts(m2, r2))
def plot_masses(ax, masses):
    uniques, counts = np.unique(masses, return_counts=True)
    ax.plot(uniques, counts)

initial_mass_factor = 1e3 / 5**-2.3 
init_5m = np.ones(1000) * 5
init_10m = np.ones(int(initial_mass_factor * 10**-2.3)) * 10
init_15m = np.ones(int(initial_mass_factor * 15**-2.3)) * 15
masses = np.concatenate((init_5m, init_10m, init_15m))



initial_radii = agn.r_g * 10**5 * np.random.uniform(0, 1, size=len(masses))**(1/3)

def num_mergers(guess, masses, initial_radii):
    mass = masses.copy()
    rad = initial_radii.copy()
    new_masses = []
    i = 0
    while i < len(mass) - 1:
        j = 1
        while i + j < len(mass) - 1:
            nm1 = len(mass)
            # nm1 = 1
            p = prob_encounter(guess, mass[i], mass[i + j], nm1, nm1, rad[i], rad[i + j])
            if p > 1 or (p > np.random.uniform(0, 1)):
                new_masses.append(mass[i] + mass[i + j])
                mass = np.delete(mass, [i, i + j])
                break
            j += 1
        i += 1
    # print(0.1 * len(masses), len(new_masses))
    # print(0.1 * len(0.1 * masses) - len(new_masses))
    return 0.1 * len(0.1 * masses) - len(new_masses)
def find_factor():
    from scipy.optimize import brentq
    factor = brentq(num_mergers, 0, 1, args=(init_5m, initial_radii[:len(init_5m)]))
    return factor
factor = find_factor()
# factor = 0.0002

fig, ax = plt.subplots()
plot_masses(ax, masses)
ax.set(yscale='log', xscale='log', xlabel='Mass ($M_\odot$)', ylabel='$N_{BH}$', ylim=[1, 1050])

for step in range(3):
    new_masses = []
    new_radii = []
    i = 0
    indices = np.random.permutation(len(masses))
    masses = masses[indices]
    initial_radii = initial_radii[indices]
    for k, radius in enumerate(initial_radii):
        ts = agn.migration_ts(masses[k], radius)
        initial_radii[k] *= abs((ts - 0.1) / ts)
    while i < len(masses) - 1:
        j = 1
        merge = 0
        while i + j < len(masses) - 1:
            nm1 = np.count_nonzero(masses == masses[i])
            nm2 = np.count_nonzero(masses == masses[i + j])
            p = prob_encounter(factor, masses[i], masses[i + j], nm1, nm2, initial_radii[i], initial_radii[i + j])
            if p > 1 or (p > np.random.uniform(0, 1)):
                new_masses.append(masses[i] + masses[i + j])
                new_radii.append(min(initial_radii[i], initial_radii[i + j]))
                masses = np.delete(masses, [i, i + j])
                initial_radii = np.delete(initial_radii, [i, i + j])
                merge = 1
                break
            j += 1
        i = i if merge else i + 1
    masses = np.concatenate((masses, new_masses))
    initial_radii = np.concatenate((initial_radii, new_radii))
    plot_masses(ax, masses)
        