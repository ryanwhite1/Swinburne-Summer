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
import scipy.interpolate as interp

G_pc = 4.3e-3   # pc.(km/s)^2/M_odot
c = 3e8     # m/s
M_odot = 1.98e30 # kg
pc_to_m = 3.086e16  # m

sigma_data_r = [0.5, 1, 1.3, 1.5, 1.7, 2, 2.6, 3, 3.5, 4, 5.5, 7]
sigma_data = [3.6, 4.1, 4.5, 4.9, 5.1, 4.8, 5.9, 5.9, 5, 4, 2, 0]
log_sigma_spline = interp.CubicSpline(sigma_data_r, sigma_data, extrapolate=True)

temp_data_r = [0.5, 1, 1.7, 2, 2.5, 3, 4, 5, 6, 7]
temp_data = [5.95, 5.75, 5.4, 5.1, 5, 4.75, 4, 3.2, 2.5, 1.8]
log_temp_spline = interp.CubicSpline(temp_data_r, temp_data, extrapolate=True)

aratio_data_r = [0.5, 0.7, 1, 1.4, 1.7, 1.9, 2, 2.2, 2.6, 3, 3.1, 3.25, 3.5, 4, 5, 6, 7]
aratio_data = [-0.7, -0.8, -0.92, -1.3, -1.45, -1.4, -1.4, -1.7, -2.1, -2.15, -2.1, -2, -1.8, -1.6, -1.05, -0.6, -0.1]
log_aratio_spline = interp.CubicSpline(aratio_data_r, aratio_data, extrapolate=True)

opacity_data_r = [0.5, 1, 1.5, 1.7, 2, 2.5, 3, 3.5, 4, 4.1, 4.2, 4.5, 5, 5.5, 6, 7]
opacity_data = [-0.4, -0.4, -0.4, -0.4, 0, -0.15, 0, -0.3, -0.4, -0.38, -0.4, -2, -3.1, -3.12, -3.15, -3.15]
log_opacity_spline = interp.CubicSpline(opacity_data_r, opacity_data, extrapolate=True)

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
        self.r_g = G_pc * self.mass / 9e10 # GM/c^2 in units of pc
        self.r_s = self.r_g * 2
        
        
        self.lenscale = lenscale    # pc
        self.lenscale_m = lenscale * pc_to_m   # m
        self.lenscale_rs = self.lenscale / self.r_s
        self.nondim_rs = self.r_s / self.lenscale
        self.timescale = np.sqrt(self.lenscale_m**3 / (G_pc * self.mass * pc_to_m * 1000**2))   # s
        self.velscale = np.sqrt(G_pc * self.mass / self.lenscale) * 1000  # m/s
        self.massscale = self.mass * M_odot    # kg
    
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
        # radial_vec2 = np.array([x, y, z]) / radius
        theta_vec = np.array([-y, x, 0]) / radius
        # theta = np.arctan2(y, x) + np.pi
        # phi = np.sign(y) * np.arccos(x / np.sqrt(x**2 + y**2))
        # radial_vec = np.array([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)])
        # print(radial_vec)
        # print(radial_vec - radial_vec2)
        # theta_vec = np.array([-np.sin(theta), np.cos(theta), 0])
        nondim_mig = self.mig_force(mass, radius) / mass * theta_vec
        nondim_damp = self.damp_force(mass, position, vel)
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
        stefboltz = 5.67 * 10**-8 * self.lenscale_m**2 * self.timescale     # units of J/K^4
        q = mass    # mass ratio of the migrator to the SMBH - assume SMBH has a mass of 1
        gamma = 5/3     # adiabatic index
        c_v = 14304 * self.massscale    # specific heat capacity of hydrogen H2 gas, units are J/K
        # logr = np.log10(radius * self.lenscale_rs)  # find log(radius) where radius is in units of SMBH Schwarzschild radii
        logr = np.log10(radius / self.nondim_rs)
        Sigma = self.disk_surfdens(logr)        # surface density
        rotvel = self.disk_rotvel(logr)
        asp_ratio = self.disk_aspectratio(logr)
        tau = self.disk_opacity(logr) * Sigma / 2  # tau = kappa * Sigma / 2
        # tau = self.disk_optdepth(logr)
        tau_eff = 3 * tau / 8 + np.sqrt(3) / 4 + 1 / (4 * tau)      # effective optical depth
        # print(tau, tau_eff)
        
        # alpha = - d(ln Sigma)/d(ln r)
        # if logr <= 3:
        #     alpha = - 1 * np.log(10)
        # else:
        #     alpha = 5/3 * np.log(10)
        # beta = - d(ln T)/d(ln r)
        # if logr <= 2.8:
        #     beta = 1 / 2.3 * np.log(10)
        # else:
        #     beta = 5/6 * np.log(10)
        alpha = -log_sigma_spline.derivative()(logr)
        beta = -log_temp_spline.derivative()(logr)
        xi = beta - (gamma - 1) * alpha
        
        Theta = (c_v * Sigma * rotvel * tau_eff) / (12 * np.pi * stefboltz * self.disk_temp(logr)**3)
        Gamma_0 = (q * radius / asp_ratio)**2 * Sigma * radius**4 * rotvel**2
        Gamma_iso =  -0.85 - alpha - 0.9 * beta
        Gamma_ad = (-0.85 - alpha - 1.7 * beta + 7.9 * xi / gamma) / gamma
        
        Gamma = Gamma_0 * (Gamma_ad * Theta**2 + Gamma_iso) / (Theta + 1)**2
        
        # print(Gamma_0, Gamma_ad, Gamma_iso, Gamma)
        return Gamma / radius   # this is the *force*... divide by mass later!
    
    def damp_force(self, mass, position, vel):
        ''' Eccentricity damping force as in Cresswell and Nelson (2007)
        Parameters
        ----------
        mass : float
        position : 1x3 np.array
        vel : 1x3 np.array
        '''
        radius = np.linalg.norm(position)
        # logr = np.log10(radius * self.lenscale_rs)  # log(radius) in schwarzschild radii    
        logr = np.log10(radius / self.nondim_rs)
        a = semi_major_axis(position, vel)
        h = self.disk_aspectratio(logr)
        # velocity = self.disk_angularvel(position, vel)
        velocity = self.disk_rotvel(logr)
        smbhmass = 1
        tdamp = (smbhmass**2 * h**4) / (mass * self.disk_surfdens(logr) * a**2 * velocity)
        e = eccentricity(position, vel)        # https://astronomy.stackexchange.com/questions/29005/calculation-of-eccentricity-of-orbit-from-velocity-and-radius
        # print(e)
        eps = e / h
        t_e = (tdamp / 0.78) * (1 - 0.14 * eps**2 + 0.06 * eps**3)
        # print(t_e)
        f_damp = -2 * np.dot(vel, position) * position / (radius**2 * t_e)
        return f_damp
    
    
    def disk_temp(self, logr):
        ''' Returns AGN disk temperature with units of K
        Parameters
        ----------
        logr : float
            log10(radius) where the radius is in units of schwarzschild radii
        '''
        # if logr <= 2.8:
        #     return 10**(-1/2.3 * logr + 6.217)
        # else:
        #     return 10**(-5/6 * logr + 22/3)
        log_temp = log_temp_spline(logr) 
        return 10**log_temp
        
    def disk_surfdens(self, logr):
        '''AGN Disk surface density, commonly with symbol Sigma. +1 in the power to go from g/cm^2 to kg/m^2
        Returns units of kg.m^-2, then non-dimensionalised
        Parameters
        ----------
        logr : float
            log10(radius) where the radius is in units of schwarzschild radii
        '''
        # if logr <= 3:
        #     val = 10**(logr + 3 + 1)
        # else:
        #     val = 10**(-5/3 * logr + 11 + 1)
        # return val * self.lenscale_m**2 / self.massscale
        log_sigma_cgs = log_sigma_spline(logr) 
        sigma = 10**(log_sigma_cgs + 1)
        return sigma * self.lenscale_m**2 / self.massscale
    
    def disk_rotvel(self, logr):
        '''Returns units of m/s, then non-dimensionalised
        Parameters
        ----------
        logr : float
            log10(radius) where the radius is in units of schwarzschild radii
        '''
        # v = np.sqrt(4.3 * 10**-3 * self.mass / (10**logr * 2 * self.r_g)) * 1000
        # # print('rotvel=', v)
        # return v / self.velscale
        v = np.sqrt(4.3 * 10**-3 * self.mass / (10**logr * 2 * self.r_g)**3) * 1000
        # print('rotvel=', v)
        return v * self.lenscale / self.velscale
    
    def disk_angularvel(self, position, velocity):
        ''' Angular velocity given position and velocity: omega = (||r x v|| / ||r||^2)
        '''
        return np.linalg.norm(np.cross(position, velocity)) / np.linalg.norm(position)**2
    
    def disk_aspectratio(self, logr):
        '''Aspect ratio of scale height: h = H / r. Unitless.
        Parameters
        ----------
        logr : float
            log10(radius) where the radius is in units of schwarzschild radii
        '''
        # if logr <= 3:
        #     val = 10**(-2/3 * logr - 1/3)
        # else:
        #     val = 10**(0.5 * logr - 3.8)
        # return val
        log_aratio = log_aratio_spline(logr) 
        return 10**(log_aratio)
        
    def disk_opacity(self, logr):
        '''Commonly with symbol kappa. -1 in the power to convert from cm^2/g to m^2/kg. 
        Units are in m^2/kg, then non-dimensionalised
        Parameters
        ----------
        logr : float
            log10(radius) where the radius is in units of schwarzschild radii
        '''
        # if logr <= 4.3:
        #     val = 10**(-0.5 - 1)
        # elif 4.3 < logr < 4.8:
        #     val = 10**(-5.6 * logr + 23.58 - 1)
        # else:
        #     val =  10**(-3.3 - 1)
        # return val * self.massscale / self.lenscale_m**2
        log_opacity_cgs = log_opacity_spline(logr) 
        kappa = 10**(log_opacity_cgs - 1)
        return kappa * self.massscale / self.lenscale_m**2
        
    def disk_optdepth(self, logr):
        '''Unitless?
        Parameters
        ----------
        logr : float
            log10(radius) where the radius is in units of schwarzschild radii'''
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
    Parameters
    ----------
    captured : list
        List of *indices* of already-merged BHs. Any index in this list will mean that a BH with that index wont have forces
        calculated and wont exert forces
    forces : bool
        If true, calculates forces from the AGN disk
    
    '''
    check_inds = np.setdiff1d(np.arange(len(masses)), captured)     # find which BH indices are not in the captured list
    
    # first a half-step kick
    vel[:] = vel + 0.5 * dt * accel # note that you must slice arrays to modify them in place in the function!
    # then full-step drift
    pos[:] = pos + dt * vel
    # then recompute accelerations
    accel[:] = ptg.Accel(pos, masses, softening, parallel=True)
    if forces:
        for i in check_inds:
            # b = accel[i]
            accel[i] += agn.get_forces(masses[i], pos[i], vel[i])
            # print(b)
    # then another half-step kick
    vel[:] = vel + 0.5 * dt * accel
    
    
    
    # now check for captures
    
    for i, primary in enumerate(check_inds):    # iterate over all non-merged BHs
        r1 = np.linalg.norm(pos[primary])
        m1 = masses[primary]
        for j, secondary in enumerate(check_inds[i + 1:]):  # now iterate over every other BH 
            r2 = np.linalg.norm(pos[secondary])
            m2 = masses[secondary]
            R_mH = np.cbrt((m1 + m2) / (3 * masses[0])) * (r1 + r2) / 2     # calculate mutual hill radius
            # print(pos[primary] - pos[secondary])
            dist = np.linalg.norm(pos[primary] - pos[secondary])    # calculate the distance between the two BHs at this timestep
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
                print(check_inds)
                captured.append(primary)    # add the primary to the captured list
                captured[:] = captured      # modify in place to update the list outside of this function
                masses[secondary] += masses[primary]; masses[secondary] *= 0.95   # assume merged mass is 95% of the sum of the original masses
                pos[secondary] = (m1 * pos[primary] + m2 * pos[secondary]) / (m1 + m2)  # set the position of the new BH to be the mass-weighted average
                vel[secondary] = (m1 * vel[primary] + m2 * vel[secondary]) / (m1 + m2)  # set the velocity of the new BH to be the mass-weighted average
                masses[primary] = 0 # set the mass of the 'other' BH to 0 so that it doesnt affect the rest of the sim
                break       # we want to break because we dont want to merge the same BH more than once in 1 timestep
            
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
    captured = [0]      # list of captured BHs to set them at the origin and not calculate migration, etc. Start with the SMBH 
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
    np.random.seed(seed) # seed the RNG for reproducibility
    
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
    yprop = -np.cos(angle)
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
    softening = np.repeat(0.1, N + 1) if N > 4e3 else np.repeat(0.0005, N + 1)
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
    a = radius / (2 - radius * vel_mag**2)    # https://physics.stackexchange.com/questions/295431/how-can-i-calculate-the-semi-major-axis-from-velocity-position-and-pull
    # also found by rearranging the expression v^2 = GM(2/r - 1/a)
    return a

def eccentricity(position, velocity):
    '''https://en.wikipedia.org/wiki/Eccentricity_vector
    '''
    radius = np.linalg.norm(position)
    e = np.linalg.norm(np.cross(velocity, np.cross(position, velocity)) - position / radius)
    return e
    
    

    

    
    
    
    