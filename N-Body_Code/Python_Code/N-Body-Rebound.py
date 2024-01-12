# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 12:08:43 2023

@author: ryanw
"""

import rebound
import numpy as np
import matplotlib.pyplot as plt
import time
from CommonTools import *


def rebound_ICs(sim, pos, masses, vel):
    n = len(masses)
    
    for i in range(n):
        x, y, z = pos[i]
        vx, vy, vz = vel[i]
        sim.add(m=masses[i], x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, hash=i)
        
def disk_forces(reb_sim):
    for particle in sim.particles[1:]:
        pos = np.array([particle.x, particle.y, particle.z])
        vel = np.array([particle.vx, particle.vy, particle.vz])
        accel = agn.get_forces(particle.m, pos, vel)
        # print(particle.hash)
        # print(accel)
        # print([particle.ax, particle.ay, particle.az])
        particle.ax += accel[0]; particle.ay += accel[1]; particle.az += accel[2]
        # print([particle.ax, particle.ay, particle.az])
    
    for i, primary in enumerate(sim.particles):    # iterate over all non-merged BHs
        if i == 0: continue
        primary_pos = np.array([primary.x, primary.y, primary.z])
        r1 = np.linalg.norm(primary_pos)
        m1 = primary.m
        for j, secondary in enumerate(sim.particles[i + 1:]):  # now iterate over every other BH 
            secondary_pos = np.array([secondary.x, secondary.y, secondary.z])
            r2 = np.linalg.norm(secondary_pos)
            m2 = secondary.m
            R_mH = np.cbrt((m1 + m2) / (3 * sim.particles[0].m)) * (r1 + r2) / 2     # calculate mutual hill radius
            # print(pos[primary] - pos[secondary])
            dist = np.linalg.norm(primary_pos - secondary_pos)    # calculate the distance between the two BHs at this timestep
            # print(R_mH)
            if dist < R_mH: # check if they within their mutual hill radius
                ### Below commented out lines check for relative kinetic energy vs binding energy in capture
                # reduced_mass = 1 / (1 / m1 + 1 / m2)
                # rel_kin_energy = 0.5 * reduced_mass * np.linalg.norm(vel[primary] - vel[secondary])**2
                # binding_energy = m1 * m2 / (2 * R_mH)
                # print(rel_kin_energy, binding_energy)
                # if rel_kin_energy < binding_energy:
                #     print("capture!", R_mH, dist)
                #     print(check_inds)
                #     captured.append(primary)    # add the primary to the captured list
                #     captured[:] = captured      # modify in place to update the list outside of this function
                #     masses[secondary] += masses[primary]    # assume no mass loss in the merger
                #     masses[primary] = 0 # set the mass of the 'other' BH to 0 so that it doesnt affect the rest of the sim
                #     break       # we want to break because we dont want to merge the same BH more than once in 1 timestep
            
                print("capture!", R_mH, dist)
                # print(check_inds)
                # captured.append(primary)    # add the primary to the captured list
                # captured[:] = captured      # modify in place to update the list outside of this function
                # masses[secondary] += masses[primary]; masses[secondary] *= 0.95   # assume merged mass is 95% of the sum of the original masses
                # pos[secondary] = (m1 * pos[primary] + m2 * pos[secondary]) / (m1 + m2)  # set the position of the new BH to be the mass-weighted average
                # vel[secondary] = (m1 * vel[primary] + m2 * vel[secondary]) / (m1 + m2)  # set the velocity of the new BH to be the mass-weighted average
                # masses[primary] = 0 # set the mass of the 'other' BH to 0 so that it doesnt affect the rest of the sim
                
                primary.m = 0.95 * (m1 + m2)
                pos = (m1 * primary_pos + m2 * secondary_pos) / (m1 + m2)
                primary.x = pos[0]; primary.y = pos[1]; primary.z = pos[2]
                primary_vel = np.array([primary.vx, primary.vy, primary.vz]); secondary_vel = np.array([secondary.vx, secondary.vy, secondary.vz])
                vel = (m1 * primary_vel + m2 * secondary_vel) / (m1 + m2)
                primary.vx = vel[0]; primary.vy = vel[1]; primary.vz = vel[2]
                sim.remove(i + j)
                break       # we want to break because we dont want to merge the same BH more than once in 1 timestep

def run_sim(sim, Tmax, dt):
    # times = np.arange(0, Tmax, dt)
    # for i, time in enumerate(times):
        
    #     sim.step()
    sim.integrate(Tmax)
        
        

Tmax = 20000
dt = 0.05
nt = int((Tmax - 0) / dt) + 1
NsBH = 100
NsBHMasses = np.ones(NsBH) * 10
SMBHMass = 1e8
r_s = 2 * 4.3e-3 * SMBHMass / 9e10  # 2GM / c^2     units of pc
Nr_s = 1e3      # number of schwarzschild radii to initialise the sim with respect to
lenscale = Nr_s * r_s
agn = AGNDisk(SMBHMass, lenscale)

pos, masses, vel, softening = AGNBHICs(NsBHMasses, SMBHMass)

sim = rebound.Simulation()
sim.integrator = "ias15"
sim.additional_forces = disk_forces
sim.force_is_velocity_dependent = 1
sim.softening = 0.1

rebound_ICs(sim, pos, masses, vel)

# sim.integrate(10000)
sim.start_server(port=1234)
run_sim(sim, 10000, dt)
op = rebound.OrbitPlot(sim)
