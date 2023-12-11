# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 07:23:22 2023

@author: ryanw
"""

import pytreegrav as ptg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import time

from CommonTools import *

Tmax = 20000
dt = 0.1
nt = int((Tmax - 0) / dt) + 1
NsBH = 10
NsBHMasses = np.ones(NsBH) * 10
SMBHMass = 1e8
r_s = 2 * 4.3e-3 * SMBHMass / 9e10  # 2GM / c^2     units of pc
Nr_s = 1e3      # number of schwarzschild radii to initialise the sim with respect to
lenscale = Nr_s * r_s
agn = AGNDisk(SMBHMass, lenscale)

## --- Disk Model --- ###
n = 2000
radius = np.logspace(0.5, 7, n) * agn.nondim_rs
logradius = np.log10(radius / agn.nondim_rs)

fig, axes = plt.subplots(nrows=5, sharex=True, figsize=(7, 10))
temps, sigma, h, tau, kappa = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
for i in range(n):
    logr = logradius[i]
    temps[i] = agn.disk_temp(logr)
    sigma[i] = agn.disk_surfdens(logr)
    h[i] = agn.disk_aspectratio(logr)
    tau[i] = agn.disk_optdepth(logr)
    kappa[i] = agn.disk_opacity(logr)
sigma = sigma * agn.massscale / agn.lenscale_m**2
kappa = kappa * agn.lenscale_m**2 / agn.massscale
axes[0].plot(radius / agn.nondim_rs, temps)
axes[1].plot(radius / agn.nondim_rs, sigma)
axes[2].plot(radius / agn.nondim_rs, h)
axes[3].plot(radius / agn.nondim_rs, kappa)
axes[4].plot(radius / agn.nondim_rs, tau, label='approximated')
axes[4].plot(radius / agn.nondim_rs,  kappa * sigma / 2, label='calculated')

for i, ax in enumerate(axes):
    ax.set(xscale='log', yscale='log')
axes[0].set(ylabel='Temperature (K)')
axes[1].set(ylabel='Surface Density (kg/m$^2$)')
axes[2].set(ylabel='Aspect Ratio (H/r)')
axes[3].set(ylabel='Opacity (m$^2$/kg)')
axes[4].set(ylabel='Optical Depth', xlabel='log(r/R_s)')
axes[4].legend()
fig.savefig('Disk Model.png', dpi=400, bbox_inches='tight')


### --- Torques Calculation --- ###



radii = np.logspace(0.6, 5, 1000) * agn.nondim_rs
torques = np.zeros(len(radii))
for i, radius in enumerate(radii):
    torques[i] = agn.mig_force(1e-7, radius)

torques *= agn.massscale * agn.lenscale_m**2 * 1e7 / (agn.timescale**2 * 1e49)
pos_vals = torques > 0 
neg_vals = torques <= 0 
fig, ax = plt.subplots()

ax.scatter(radii[pos_vals] / agn.nondim_rs, torques[pos_vals], s=5, label='negative')
ax.scatter(radii[neg_vals] / agn.nondim_rs, -torques[neg_vals], c='r', s=5, label='positive')
ax.set(xscale='log', yscale='log', xlabel="Log(R/R$_s$)", ylabel='abs($\Gamma$) / 1e49')
ax.legend()
fig.savefig('Torque Model.png', dpi=400, bbox_inches='tight')