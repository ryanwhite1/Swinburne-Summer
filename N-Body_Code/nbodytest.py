# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 18:57:57 2023

@author: ryanw
"""

import pytreegrav as ptg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import time

from CommonTools import *


Tmax = 10000
dt = 0.1
nt = int((Tmax - 0) / dt) + 1
NsBH = 10
NsBHMasses = np.ones(NsBH) * 10
SMBHMass = 4e6
r_s = 2 * 4.3 * 10**-3 * SMBHMass / 9e10
Nr_s = 1e4
lenscale = Nr_s * r_s
agn = AGNDisk(SMBHMass, lenscale)

pos, masses, vel, softening = AGNBHICs(NsBHMasses, SMBHMass)
positions, velocities = perform_sim(Tmax, dt, pos, masses, vel, softening, agn, forces=True)

times = np.linspace(0, Tmax, nt)
real_times = time_convert(times, SMBHMass + sum(NsBHMasses), lenscale)

# animate_sim(positions, 'mig_test', 20, every=10, times=[True, real_times])


fig, ax = plt.subplots()
step = 5
for i in range(1, 11):
    radii = np.array([np.linalg.norm(positions[i, :, j]) for j in range(0, nt, step)])
    vel_mags = np.array([np.linalg.norm(velocities[i, :, j]) for j in range(0, nt, step)])
    semi_majors = - radii / (radii * vel_mags**2 - 2)
    ax.plot(real_times[::step], semi_majors * Nr_s)
ax.set(yscale='log', xlabel="Time (Myr)", ylabel="Semi-Major Axis ($R_s$)")
# fig.savefig('NBodyTest.png', dpi=400, bbox_inches='tight')
# fig.savefig('MigrationTest.png', dpi=400, bbox_inches='tight')
fig.savefig('CaptureTest.png', dpi=400, bbox_inches='tight')