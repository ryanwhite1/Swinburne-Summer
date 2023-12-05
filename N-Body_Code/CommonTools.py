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

G_pc = 4.3e-3   # pc.(km/s)^2/M_odot
c = 3e8     # m/s

class AGNDisk(object):
    def __init__(self, smbhmass, lenscale):
        '''
        Parameters
        ----------
        smbhmass : float
            Mass of the SMBH in units of solar masses.
        lenscale : float
            Characteristic length scale of the simulation (units of pc)
        '''
        self.mass = smbhmass    # sol masses
        # self.mass_ratio = smbhmass / 1e8
        self.r_g = 4.3e-3 * self.mass / 9e10 # GM/c^2 in units of pc
        
        self.lenscale = lenscale    # pc
        self.lenscale_m = lenscale * 3.086e16   # m
        self.lenscale_rs = self.lenscale / (2 * self.r_g)
        self.timescale = np.sqrt(self.lenscale_m**3 / (4.3e-3 * self.mass * 3.086e16 * 1000**2))   # s
        self.velscale = np.sqrt(4.3e-3 * self.mass / self.lenscale) * 1000  # m/s
        self.massscale = self.mass * 1.98e30    # kg
    
    
    # def migration_ts(self, mass, radius):
    #     '''McKernan et al (2018)'''
    #     surf_dens_ratio = 10**5 / self.disk_surfdens(radius) # (sigma / 10^5)^-1
    #     mass_ratio = 5 / mass # (mass / 5M_odot)^-1
    #     scale_height_ratio = (self.disk_aspectratio(radius) / 0.02)**2 # assume center of mass radius is just orbital radius, i.e. M_smbh >> M_bh
    #     orbit_ratio = (radius / (10**4 * self.r_g))**(-1/2)
    #     ts = 38 * orbit_ratio * mass_ratio * scale_height_ratio * surf_dens_ratio * self.mass_ratio**(3/2)
    #     return ts
    
    def get_forces(self, mass, position, vel):
        '''
        Parameters
        ----------
        mass : float
            Mass of the migrating object
        position : (1x3) np.array
            Position (x, y, z) in units of Schwarzschild radii from the origin SMBH
        # radius : float
        #     Radius from the SMBH in units of its Schwarzschild radius.
        '''
        x, y, z = position
        radius = np.linalg.norm(position)
        radial_vec = np.array([x, y, z]) / radius
        theta_vec = np.array([y, -x, 0]) / radius
        nondim_mig = self.mig_force(mass, radius) * theta_vec
        nondim_damp = self.damp_force(mass, position, vel) * radial_vec
        # print(nondim_mig, nondim_damp)
        
        return nondim_mig + nondim_damp
    
    def mig_force(self, mass, radius):
        '''
        Parameters
        ----------
        mass : float
            Mass of the migrating object in N-body units.
        radius : float
            Radius from the SMBH in N-Body units.
        '''
        stefboltz = 5.67 * 10**-8 * self.lenscale_m**2 * self.timescale
        # q = mass / self.mass    # mass ratio of the migrator to the SMBH
        q = mass
        # pc_to_m = 3.086 * 1e16
        gamma = 5/3     # adiabatic index
        c_v = 14304 * self.massscale    # specific heat capacity of hydrogen H2 gas
        
        logr = np.log10(radius * self.lenscale_rs)  # find log(radius) where radius is in units of SMBH Schwarzschild radii
        # radius_m = radius * self.lenscale_m
        Sigma = self.disk_surfdens(logr)
        rotvel = self.disk_rotvel(logr)
        asp_ratio = self.disk_aspectratio(logr)
        tau = self.disk_opacity(logr) * Sigma / 2  # tau = kappa * Sigma / 2
        # tau = self.disk_optdepth(logr)
        tau_eff = 3 * tau / 8 + np.sqrt(3) / 4 + 1 / (4 * tau)
        # print(tau, tau_eff)
        
        # alpha = - d(ln Sigma)/d(ln r)
        if logr <= 3:
            alpha = - np.log(10)
        else:
            alpha = np.log(10) * 5/3 
        # beta = - d(ln T)/d(ln r)
        if logr <= 2.8:
            beta = np.log(10) / 2.3 
        else:
            beta = np.log(10) * 5/6
        xi = beta - (gamma - 1) * alpha
        
        Theta = (c_v * Sigma * rotvel * tau_eff) / (12 * np.pi * stefboltz * self.disk_temp(logr)**3)
        # print(Theta)
        Gamma_0 = (q / asp_ratio)**2 * Sigma * radius**4 * rotvel**2
        # print(Gamma_0, (q / asp_ratio), Sigma, radius**4, rotvel**2)
        Gamma_iso =  -0.85 - alpha - 0.9 * beta
        Gamma_ad = (-0.85 - alpha - 1.7 * beta + 7.9 * xi / gamma) / gamma
        
        Gamma = Gamma_0 * (Gamma_ad * Theta**2 + Gamma_iso) / (Theta + 1)**2
        # print(Gamma_0, Gamma_ad, Gamma_iso, Gamma)
        return Gamma / (radius * mass)
    
    def damp_force(self, mass, position, vel):
        '''
        Parameters
        ----------
        radius : float
            Radius from the SMBH in N-Body units.
        '''
        radius = np.linalg.norm(position)
        logr = np.log10(radius * self.lenscale_rs)  # log(radius) in schwarzschild radii    
        a = semi_major_axis(position, vel)
        h = self.disk_aspectratio(logr)
        smbhmass = 1
        tdamp = (smbhmass**2 * h**4) / (mass * smbhmass * self.disk_surfdens(logr) * a**2 * self.disk_rotvel(logr))
        e = np.linalg.norm(np.cross(vel, np.cross(position, vel)) - position / radius)        # https://astronomy.stackexchange.com/questions/29005/calculation-of-eccentricity-of-orbit-from-velocity-and-radius
        # print(e)
        eps = e / h
        t_e = (tdamp / 0.78) * (1 - 0.14 * eps**2 + 0.06 * eps**3)
        # print(t_e)
        f_damp = -2 * np.dot(vel, position) * position / (radius**2 * t_e)
        return f_damp
    
    
    def disk_temp(self, logr):
        ''' Returns units of K'''
        if logr <= 2.8:
            return 10**(-1/2.3 * logr + 6.217)
        else:
            return 10**(-5/6 * logr + 22/3)
    def disk_surfdens(self, logr):
        '''+1 in the power to go from g/cm^2 to kg/m^2
        Returns units of kg.m^-2'''
        if logr <= 3:
            val = 10**(logr + 3 + 1)
        else:
            val = 10**(-5/3 * logr + 11 + 1)
        return val * self.lenscale_m**2 / self.massscale
    def disk_rotvel(self, logr):
        '''Returns units of m/s'''
        # v = np.sqrt(4.3 * 10**-3 * self.mass / (10**logr * 2 * self.r_g)) * 1000
        # # print('rotvel=', v)
        # return v / self.velscale
        v = np.sqrt(4.3 * 10**-3 * self.mass / (10**logr * 2 * self.r_g)**3) * 1000
        # print('rotvel=', v)
        return v * self.lenscale / self.velscale
    def disk_aspectratio(self, logr):
        '''Aspect ratio of scale height: h = H / r. Unitless.'''
        if logr <= 3:
            val = 10**(-2/3 * logr - 1/3)
        else:
            val = 10**(0.5 * logr - 3.5)
        return val
        
    def disk_opacity(self, logr):
        '''-1 in the power to convert from cm^2/g to m^2/kg. Units are in m^2/kg'''
        if logr <= 4.3:
            val = 10**(-0.5 - 1)
        elif 4.3 < logr < 4.8:
            val = 10**(-5.6 * logr + 23.58 - 1)
        else:
            val =  10**(-3.3 - 1)
        return val * self.massscale / self.lenscale_m**2
        
    def disk_optdepth(self, logr):
        '''Unitless?'''
        if logr <= 3:
            return 10**(logr + 2)
        else:
            return 10**(-2.5 * logr + 12.5)
        
        
    # def disk_toomre(self, radius):
    #     ''''''
    #     logr = np.log10(radius)
    #     if logr >= -2:
    #         return 1 
    #     else:
    #         return 10**(-5 * logr - 10)
    
    

def leapfrog_kdk_timestep(dt, pos, masses, softening, vel, accel, agn, captured, forces=True):
    '''
    '''
    check_inds = np.setdiff1d(np.arange(len(masses)), captured)
    
    # first a half-step kick
    vel[:] = vel + 0.5 * dt * accel # note that you must slice arrays to modify them in place in the function!
    # then full-step drift
    pos[:] = pos + dt * vel
    # then recompute accelerations
    accel[:] = ptg.Accel(pos, masses, softening, parallel=True)
    if forces:
        # accel[1:] += np.array([agn.get_forces(masses[i], pos[i], vel[i]) for i in range(1, len(pos[:]))])
        for i in check_inds:
            # b = accel[i]
            accel[i] += agn.get_forces(masses[i], pos[i], vel[i])
            # print(b - accel[i])
    # then another half-step kick
    vel[:] = vel + 0.5 * dt * accel
    
    
    
    # now check for captures
    
    for i, primary in enumerate(check_inds):
        r1 = np.linalg.norm(pos[primary])
        m1 = masses[primary]
        for j, secondary in enumerate(check_inds[i + 1:]):
            r2 = np.linalg.norm(pos[secondary])
            m2 = masses[secondary]
            R_mH = np.cbrt((m1 + m2) / (3 * masses[0])) * (r1 + r2) / 2 
            # print(pos[primary] - pos[secondary])
            dist = np.linalg.norm(pos[primary] - pos[secondary])
            if dist < R_mH:
                print("capture!", R_mH, dist)
                print(check_inds)
                captured.append(primary)
                captured[:] = captured
                masses[secondary] += masses[primary]
                masses[primary] = 0
                break
            
    # now set the central SMBH/captured BHs to not have changed position or velocity
    for i in captured:
        vel[i] = [0, 0, 0]; pos[i] = [0, 0, 0]; accel[i] = [0, 0, 0]
            
    
def perform_sim(Tmax, dt, pos, masses, vel, softening, agn, forces=True):
    ''' Performs an nbody sim from t=0 to t=Tmax [in steps of dt] given particles with parameters
    pos : (N x 3) ndarray
        Positions of N particles in xyz coordinates
    masses : (N x 1) array
        Mass of each particle
    vel : (N x 3) ndarray
        Velocities of N particles in xyz components
    softening : (N x 1) array
        Softening coefficient of each particle
    forces : bool
        Whether or not to model disk forces (migration, damping) in the 1d AGN disk.
    Returns
    -------
    positions : (N x 3 x nt) ndarray
        Particle position in xyz space for each of the N particles at each of the nt time steps. 
    '''
    t1 = time.time()
    N = len(pos[:, 0])
    accel = ptg.Accel(pos, masses, softening, parallel=True) # initialize acceleration

    t = 0 # initial time
    nt = int((Tmax - t) / dt) + 1


    positions = np.zeros((N, 3, nt))
    positions[:, :, 0] = pos
    velocities = np.zeros((N, 3, nt))
    velocities[:, :, 0] = vel
    
    i = 0
    captured = [0]
    while t <= Tmax: # actual simulation loop - this may take a couple minutes to run
        leapfrog_kdk_timestep(dt, pos, masses, softening, vel, accel, agn, captured, forces=forces)
        positions[:, :, i] = pos
        velocities[:, :, i] = vel
        t += dt
        i += 1
    t2 = time.time()
    print(f"Simulation complete in {round(t2 - t1, 3)}s!")
    return positions, velocities

def AGNBHICs(masses, smbhmass, seed=4080):
    ''' Initial conditions for a uniform collapse of N particles with initial velocity of vel_prop of the equilibrium velocity.
    Parameters
    ----------
    masses : np.array
        Masses of each of the initial sBHs (in M_\odot)
    smbhmass : float
        Mass (in M_\odot) of the central SMBH
    '''
    N = len(masses)
    # r_s = 4.3 * 10**-3 * smbhmass / (9 * 10**10)
    # np.random.seed(seed) # seed the RNG for reproducibility
    
    # uniformly (and randomly) distribute points in the unit disk
    theta = np.random.uniform(0, 2*np.pi, N)
    dists = np.random.uniform(0.5**3, 1, N)
    R = 1 * np.cbrt(dists)
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    z = 0
    pos = np.zeros((N, 3))
    pos[:, 0] = x; pos[:, 1] = y; pos[:, 2] = z
    # pos -= np.average(pos, axis=0) # apply small correction to put center of mass at the origin
    
    # now to give particles random velocity with a magnitude of 'vel_prop' of the equilibrium velocity
    des_vel = np.sqrt(1 / R)
    angle = np.arctan2(y, x)
    xprop = np.sin(angle)
    yprop = - np.cos(angle)
    zprop = np.zeros(N)
    mult = np.sqrt(des_vel**2 / (xprop**2 + yprop**2 + zprop**2))
    xprop *= mult; yprop *= mult
    
    vel = np.zeros_like(pos) # initialize at rest
    vel[:, 0] = xprop; vel[:, 1] = yprop; vel[:, 2] = zprop
    # vel *= np.random.choice([1, -1], size=N)
    for i in range(N):
        mag = np.sqrt(vel[i, 0]**2 + vel[i, 1]**2)
        print(des_vel[i] - mag)
    print(vel)
    # vel -= np.average(vel, axis=0) # make average velocity 0
    # insert the SMBH into the start of the array
    vel = np.insert(vel, 0, [0, 0, 0], axis=0)
    pos = np.insert(pos, 0, [0, 0, 0], axis=0)
    softening = np.repeat(0.2, N + 1) if N > 4e3 else np.repeat(0.1, N + 1)
    masses = masses / (sum(masses) + smbhmass)
    masses = np.insert(masses, 0, smbhmass / (sum(masses) + smbhmass))
    
    return pos, masses, vel, softening


    
    
def animate_sim(positions, filename, length, colours=[], every=1, times=[False]):
    ''' Animates the positions of N points in 3D space for nt timesteps against a black background, and saves it too!
    Parameters
    ----------
    positions : (N x 3 x nt) ndarray
        Particle position in xyz space for each of the N particles at each of the nt time steps. 
    filename : str
        The desired filename, to be saved as 'filename.gif'
    length : float
        Desired length (in seconds) of the gif
    colours : Nx1 list/array (optional)
        The (order dependent) colours to plot each of the N data points
    every : int
        Will plot every n frames in the animation. every=1 corresponds to plotting each frame, every=2 each second frame, etc.
    times : list
        Either 1 element list (containing just False) if we don't want a little timer in the top right, or, 
        a 2 element list (the first being True) with the second element containing an nt x 1 array of times. 
    '''
    fig = plt.figure(figsize=(12, 12), frameon=False)   # we want no frame so that it's a clean black animation
    ax = fig.add_subplot(projection='3d')   # 3d axes
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None) # removes (most of the) blank border around the plot
    # now do an initial scatter so we can get the axis limits to use through the animation
    ax.scatter(positions[:, 0, 0], positions[:, 1, 0], positions[:, 2, 0], s=1, marker='.')
    xmin, xmax = min(positions[:, 0, 0]), max(positions[:, 0, 0])
    ymin, ymax = min(positions[:, 1, 0]), max(positions[:, 1, 0])
    zmin, zmax = min(positions[:, 2, 0]), max(positions[:, 2, 0])
    limmin, limmax = min([xmin, ymin, zmin]), max([xmax, ymax, zmax])   # get the minimum and maximums for the axis limits
    
    # now calculate some parameters for the animation frames and timing
    nt = len(positions[0, 0, :]) # number of timesteps
    frames = np.arange(0, nt, every)    # iterable for the animation function. Chooses which frames (indices) to animate.
    fps = len(frames) // length  # fps for the final animation
    
    ax.set_facecolor('k')   # black background, since space is blach duh
    scales = np.ones(len(positions[:, 0, 0])) * 5 
    scales[0] = 25
        
    def animate(i):
        if (i // every)%20 == 0:
            print(f"{i // every} / {len(frames)}")
            
        ax.clear()
        if len(colours) != 0:
            ax.scatter(positions[:, 0, i], positions[:, 1, i], positions[:, 2, i], s=scales, marker='.', 
                        c=colours)
        else:
            ax.scatter(positions[:, 0, i], positions[:, 1, i], positions[:, 2, i], s=scales, marker='.', c='w')
        ax.set_xlim(limmin, limmax); ax.set_ylim(limmin, limmax); ax.set_zlim(limmin, limmax)
        if times[0]:    # plot the current time in the corner if we want to!
            ax.text(0.7 * limmax, 0.9 * limmax, 0, "$T = " + str(round(times[1][i], 2)) + "$", fontsize=24, color='w')
        
        ax.grid(False)
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.elev = 90    # sets the viewing angle to be top-down
        return fig,

    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=1, cache_frame_data=False, blit=True)
    ani.save(f"{filename}.gif", writer='pillow', fps=fps)
    plt.close('all')
    
def time_convert(time, M, R):
    ''' Converts from n-body time to real time (in units of Myr).
    Parameters
    ----------
    time : float/array
        N-body times
    M : float
        Mass of the system in units of solar masses
    R : float
        Radius of the system in pc
    '''
    mass = M * 1.988e30
    radius = R * 3.086e16
    G = 6.6743e-11
    Myr_sec = 31536000000000.0
    return time * np.sqrt(radius**3 / (mass * G)) / Myr_sec

def semi_major_axis(position, velocity):
    '''N-Body Units calc'''
    radius = np.linalg.norm(position)
    vel_mag = np.linalg.norm(velocity)
    a = - radius / (radius * vel_mag**2 - 2)    # https://physics.stackexchange.com/questions/295431/how-can-i-calculate-the-semi-major-axis-from-velocity-position-and-pull
    return a

def eccentricity(position, velocity):
    '''
    '''
    radius = np.linalg.norm(position)
    e = np.linalg.norm(np.cross(velocity, np.cross(position, velocity)) - position / radius)
    return e
    
    

    

    
    
    
    