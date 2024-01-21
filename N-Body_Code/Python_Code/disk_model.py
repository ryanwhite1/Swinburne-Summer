# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 11:30:52 2024

@author: ryanw
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp

# Teff = x[0], T = x[1], taueff = x[2], tau = x[3], kappa = x[4], Sigma = x[5], cs = x[6]
# rho = x[7], h = x[8]

root = [3.5, 5., 5., 6., 0., 6., 8., -10., -2.]
G_pc = 4.3e-3
G = 6.67e-11
M_odot = 1.98e30
c = 299792458
c_cgs = c * 100
mu = 0.62
m_H = 1.6735575e-24     # cgs units
m_cgs = mu * m_H
thomson_cross_sec = 6.65246e-29

def kappa_formula(rho, T, log=False):
    rho /= (T/1e6)**3
    if T <= 166.81:
        k0, a, b = 2e-4, 0, 2
    elif 166.81 < T <= 202.677:
        k0, a, b = 2e16, 0, -7
    elif 202.677 < T <= 2286.77 * rho**(2/49):
        k0, a, b = 1e-1, 0, 0.5
    elif 2286.77 * rho**(2/49) < T <= 2029.76 * rho**(1/81):
        k0, a, b = 2e81, 1, -24
    elif 2029.76 * rho**(1/81) < T <= 1e4 * rho**(1/21):
        k0, a, b = 1e-8, 2/3, 3
    elif 1e4 * rho**(1/21) < T <= 31195.2 * rho**(4/75):
        k0, a, b = 1e-36, 1/3, 10
    elif 31195.2 * rho**(4/75) < T <= 1.79393e8 * rho**(2/5):
        k0, a, b = 1.5e20, 1, -2.5
    else:
        k0, a, b = 0.348, 0, 0
    if not log:
        kap = k0 * rho**a * T**b
        if np.isnan(kap):
            kap = 0.348
        return kap
    else:
        return np.log10(k0) + a * np.log10(rho) + b * np.log10(T)

kappa_data = np.genfromtxt('X07Y029Z001.txt')
# kappa_data = np.genfromtxt('X07Y027Z003.txt')
kappa_R = kappa_data[0, 1:]
kappa_T = kappa_data[1:, 0]
kappa_interp = interp.RegularGridInterpolator((kappa_T, kappa_R), kappa_data[1:, 1:], bounds_error=False)


def kappa_from_data(logrho, logT):
    logR = logrho - (3 * (logT - 6))
    # print(logR)
    # if logR <= kappa_R[0]:
    #     if logT <= kappa_T[0]:
    #         return kappa_data[1, 1]
    #     elif logT >= kappa_T[-1]:
    #         return kappa_data[-1, 1]
    # elif logR >= kappa_R[-1]:
    #     if logT <= kappa_T[0]:
    #         return kappa_data[1, -1]
    #     elif logT >= kappa_T[-1]:
    #         return kappa_data[-1, -1]
    return kappa_interp([logT, logR])

stef_boltz = 5.67037e-5    # cgs units
k = 1.38065e-24        # cgs

def log_system(x, r, Mdot, angvel, alpha, c_cgs, m_cgs):
    return [4 * x[0] - np.log10(3 / (8 * np.pi * stef_boltz) * Mdot * angvel**2),       # Teff
            x[1] * 4 - (x[3] + x[0] * 4),                                                 # T
            x[2] - np.log10(3 * 10.**x[3] / 8 + 0.5 + 0.25 / (10.**x[3])),              # tau_eff
            x[3] - (x[4] + x[5] - np.log10(2)),                                         # tau
            # x[4] - kappa_formula(10.**x[7], 10.**x[1]),
            x[4] - kappa_from_data(x[7], x[1]),
            x[5] - (np.log10(Mdot * angvel / (3 * np.pi * alpha)) - 2 * x[6]),
            x[8] + np.log10(r) - (x[6] - np.log10(angvel)),
            x[7] - (x[5] - np.log10(2) - (x[8] + np.log10(r))),
            2 * x[6] - np.log10(k * 10.**x[1] / m_cgs + stef_boltz * 10.**x[3] * 10.**(4 * x[0]) / (2 * c_cgs * 10.**x[7]))]         


def angvel(r, M):
    return np.sqrt(G * M * M_odot / r**3)

n = 1000

M = 1e8
alpha = 0.01
f_edd = 0.1
Mdot_edd = 4 * np.pi * G * M * M_odot * (m_H / 1e3) / (0.1 * thomson_cross_sec * c)
Mdot = f_edd * Mdot_edd
rs_m = 2 * G * M * M_odot / c**2
rs_cm = rs_m * 1e2
r_min = 6 * rs_cm

log_radii = np.logspace(1.01 * np.log10(r_min / rs_cm), 6, n)
radii = log_radii * rs_cm
temps, rho, Sigma, h, kappa, tau, Q = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
t_eff = np.zeros(n)

for i, r in enumerate(radii):
    Mdotdash = Mdot * 1e3 * (1 - np.sqrt(r_min / r))
    angvel_r = angvel(r / 100, M)
    root = opt.fsolve(log_system, root, args=(r, Mdotdash, angvel_r, alpha, c_cgs, m_cgs))
    # print(np.isclose(log_system(root, Mdotdash, angvel_r, alpha, c_cgs, m_cgs), [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    temps[i] = 10**root[1]
    Sigma[i] = 10**root[5]
    rho[i] = 10**root[7]
    h[i] = 10**root[8]
    kappa[i] = 10**root[4]
    tau[i] = 10**root[3]
    Q[i] = max(h[i] * r * angvel_r**2 / (np.pi * G*1e3 * Sigma[i]), 1.4)
    
fig, axes = plt.subplots(nrows=6, sharex=True, figsize=(7, 10))

axes[0].plot(log_radii, temps)
axes[1].plot(log_radii, rho)
axes[2].plot(log_radii, h)
axes[3].plot(log_radii, kappa)
axes[4].plot(log_radii, tau, label='approximated')
axes[4].plot(log_radii, kappa * Sigma / 2, label='calculated')
axes[5].plot(log_radii, Q)

for i, ax in enumerate(axes):
    ax.set(xscale='log', yscale='log')
axes[0].set(ylabel='Temperature (K)')
axes[1].set(ylabel='Surface Density (kg/m$^2$)')
axes[2].set(ylabel='Aspect Ratio (H/r)')
axes[3].set(ylabel='Opacity (m$^2$/kg)')
axes[4].set(ylabel='Optical Depth')
axes[5].set(ylabel='Toomre', xlabel='log(r/R_s)')
axes[4].legend()

















