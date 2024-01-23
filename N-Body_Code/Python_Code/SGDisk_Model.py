# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 11:30:52 2024

@author: ryanw
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp

# Teff = x[0], T = x[1], tau = x[2], kappa = x[3], Sigma = x[4], cs = x[5], rho = x[6], h = x[7], Q = x[8], beta = x[9], prad = x[10], pgas = x[11]

root = np.array([6, 6., 5., -1., 6., 8., -10., -2., 4., -3., 10., 4.])
# root = 10**np.array([5, 6., 5., -1., 6., 6., -9., -1., 4., 0., 3., 4.])

G_pc = 4.3e-3
G = 6.67e-11
G_cgs = G * 1e3
M_odot = 1.98e30
M_odot_cgs = M_odot * 1e2
c = 299792458.
c_cgs = c * 100
mu = 0.62
m_H = 1.6735575e-24     # cgs units
m_cgs = mu * m_H
thomson_cross_sec = 6.65246e-29
stef_boltz = 5.67037e-5    # cgs units
k = 1.38065e-16        # cgs

def kappa_formula(rho, T, log=False):
    # rho /= (T/1e6)**3
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
# kappa_data = np.genfromtxt('A96X070Y027Z003.txt')
kappa_R = kappa_data[0, 1:]
kappa_T = kappa_data[1:, 0]
kappa_interp = interp.RegularGridInterpolator((kappa_T, kappa_R), kappa_data[1:, 1:], bounds_error=False)
kappa_data = kappa_data[1:, 1:]

lowtemp_kappa_data = np.genfromtxt('lowtempX07Y028Z002.txt')[::-1, :]
low_kappa_R = lowtemp_kappa_data[-1, 1:]
low_kappa_T = lowtemp_kappa_data[:-1, 0]
low_kappa_interp = interp.RegularGridInterpolator((low_kappa_T, low_kappa_R), lowtemp_kappa_data[1:, :-1], bounds_error=False)
lowtemp_kappa_data = lowtemp_kappa_data[:-1, 1:]

kappa_T_overlap = np.intersect1d(low_kappa_T, kappa_T)

# def kappa_from_data(logrho, logT):
#     if logT <= 3:
#         return -0.76
#     logR = logrho - (3. * (logT - 6.))
#     # print(logR)
#     if logR <= kappa_R[0]:
#         if logT <= kappa_T[0]:
#             return kappa_data[1, 1]
#         elif logT >= kappa_T[-1]:
#             return kappa_data[-1, 1]
#         else:
#             return kappa_interp([logT, kappa_R[0]])
#     elif logR >= kappa_R[-1]:
#         if logT <= kappa_T[0]:
#             return kappa_data[1, -1]
#         elif logT >= kappa_T[-1]:
#             return kappa_data[-1, -1]
#         else:
#             return kappa_interp([logT, kappa_R[-1]])
#     elif logT >= kappa_T[-1]:
#         return kappa_interp([kappa_T[-1], logR])
#     elif logT <= kappa_T[0]:
#         return kappa_interp([kappa_T[0], logR])
#     return kappa_interp([logT, logR])

def kappa_from_data(logrho, logT):
    if logT < 3:
        return 0.76
    logR = logrho - (3. * (logT - 6.))
    
    # if logT <= low_kappa_T[-1]:
    #     if logT >= kappa_T[0]:
    #         prop = np.interp(logT, kappa_T_overlap, np.linspace(0, 1, len(kappa_T_overlap)))
    #         return prop * kappa_interp([logT, logR]) + (1 - prop) * low_kappa_interp([logT, logR])
    #     return low_kappa_interp([logT, logR])
    # elif logT <= kappa_T[-1]:
    #     return kappa_interp([logT, logR])
    # return kappa_interp([logT, logR])
    if logT < kappa_T[0]:
        if logR <= low_kappa_R[0]:
            if logT <= low_kappa_T[0]:
                return lowtemp_kappa_data[0, 0]
            else:
                return low_kappa_interp([logT, low_kappa_R[0]])
        elif logR >= low_kappa_R[-1]:
            if logT <= low_kappa_T[0]:
                return lowtemp_kappa_data[1, -1]
            else:
                return low_kappa_interp([logT, low_kappa_R[-1]])
        elif logT <= low_kappa_T[0]:
            return low_kappa_interp([low_kappa_T[0], logR])
        return low_kappa_interp([logT, logR])
    else:
        if logR <= kappa_R[0]:
            if logT >= kappa_T[-1]:
                return kappa_data[-1, 0]
            else:
                return kappa_interp([logT, kappa_R[0]])
        elif logR >= kappa_R[-1]:
            if logT >= kappa_T[-1]:
                return kappa_data[-1, -1]
            else:
                return kappa_interp([logT, kappa_R[-1]])
        elif logT >= kappa_T[-1]:
            return kappa_interp([kappa_T[-1], logR])
        return kappa_interp([logT, logR])

b = 0
min_Q = 1.

def log_system(x, r, Mdot, angvel, alpha, c_cgs, m_cgs):
    # Teff = x[0], T = x[1], tau = x[2], kappa = x[3], Sigma = x[4], cs = x[5], rho = x[6], h = x[7], Q = x[8], beta = x[9], prad = x[10], pgas = x[11]
    eq1 = 4 * x[0] - np.log10(3 * Mdot * angvel**2 / (8 * np.pi * stef_boltz))
    eq2 = 4 * x[1] - (4 * x[0] + np.log10(3 * 10.**x[2] / 8 + 0.5 + 1 / (4 * 10.**x[2])))
    eq3 = x[2] - (x[3] + x[4] - np.log10(2))
    eq4 = b * x[9] + 2 * x[5] + x[4] - np.log10(Mdot * angvel / (3 * np.pi * alpha))
    eq5 = x[10] - (x[2] + 4 * x[0] + np.log10(stef_boltz / (2 * c_cgs)))
    eq6 = x[11] - (x[6] + x[1] + np.log10(k / m_cgs))
    eq7 = x[9] - (x[11] - np.log10(10.**x[10] + 10.**x[11]))
    eq8 = x[4] - (np.log10(2) + x[6] + x[7] + np.log10(r))
    eq9 = x[7] + np.log10(r) - (x[5] - np.log10(angvel))
    eq10 = 2 * x[5] - (np.log10(10.**x[10] + 10.**x[11]) - x[6])
    eq11 = x[3] - kappa_from_data(x[6], x[1])
    # eq11 = x[3] - kappa_formula(10.**x[6], 10.**x[1], log=True)
    eq12 = x[6] - (np.log10(angvel**2 / (2 * np.pi * G_cgs)) - max(x[8], np.log10(min_Q)))
    if x[8] <= np.log10(min_Q):
        eq12 = x[6] - np.log10(angvel**2 / (2 * np.pi * G_cgs * min_Q))
        eq1 = 4 * x[1] - (4 * x[0] + np.log10(3 * 10.**x[2] / 8 + 0.5 + 1 / (4 * 10.**x[2])))
    else:
        eq12 = x[6] - (np.log10(angvel**2 / (2 * np.pi * G_cgs)) - x[8])
    return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12]            


def angvel(r, M):
    return np.sqrt(G_cgs * M * M_odot_cgs / r**3)

n = 100

M = 1e6
alpha = 1e-2
f_edd = 0.5
Mdot_edd = 4 * np.pi * G * M * M_odot * (m_H / 1e3) / (0.1 * thomson_cross_sec * c)
Mdot_edd = 4 * np.pi * G_cgs * M * M_odot_cgs * m_H / (0.1 * thomson_cross_sec * 1e4 * c_cgs)
Mdot = f_edd * Mdot_edd
rs_m = 2 * G * M * M_odot / c**2
rs_cm = 2 * G_cgs * M * M_odot_cgs / c_cgs**2
r_min = 6 * rs_cm

log_radii = np.logspace(1.1 * np.log10(r_min / rs_cm), 6, n)
radii = log_radii * rs_cm
temps, rho, Sigma, h, kappa, tau, Q = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
t_eff = np.zeros(n)

for i, r in enumerate(radii):
    Mdotdash = Mdot * (1 - np.sqrt(r_min / r))
    angvel_r = angvel(r, M)
    root = opt.fsolve(log_system, root, args=(r, Mdotdash, angvel_r, alpha, c_cgs, m_cgs))
    temps[i] = 10**root[1]
    Sigma[i] = 10**root[4]
    rho[i] = 10**root[6]
    h[i] = 10**root[7]
    kappa[i] = 10**root[3]
    tau[i] = 10**root[2]
    Q[i] = max(10.**root[8], min_Q)
    
    
    # root = opt.fsolve(system, root, args=(r, Mdotdash, angvel_r, alpha, c_cgs, m_cgs))
    # temps[i] = root[1]
    # Sigma[i] = root[4]
    # rho[i] = root[6]
    # h[i] = root[7]
    # kappa[i] = root[3]
    # tau[i] = root[2]
    # Q[i] = max(root[8], min_Q)
    
    
fig, axes = plt.subplots(nrows=6, sharex=True, figsize=(7, 10))

axes[0].plot(log_radii, temps)
axes[1].plot(log_radii, Sigma)
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

















