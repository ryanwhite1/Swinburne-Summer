# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 11:30:52 2024

@author: ryanw
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
import warnings
warnings.filterwarnings("ignore")

# set LaTeX font for our figures
# plt.rcParams.update({"text.usetex": True})
# plt.rcParams['font.family'] = 'serif'
# plt.rcParams['mathtext.fontset'] = 'cm'

G = 6.67e-11
G_pc = 4.3e-3
G_cgs = G * 1e3
M_odot = 1.98e30
M_odot_cgs = M_odot * 1e3
c = 299792458.
c_cgs = c * 100
mu = 0.62               # average molecular weight ?
m_H = 1.6735575e-24     # hydrogen mass, cgs units
m_cgs = mu * m_H
thomson_cross_sec = 6.65246e-29     # SI units
thomson_cgs = thomson_cross_sec * 1e4
stef_boltz = 5.67037e-5    # cgs units
k = 1.38065e-16        # cgs


def kappa_from_formula(rho, T, log=False, prescription='derdzinski'):
    if prescription == 'derdzinski':
        # https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.4522D/abstract
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
    elif prescription == 'metzger':
        # https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.3200M/abstract
        Z = 0.02
        X = 0.7
        Y = 1 - (X + Z)
        ke = 0.2 * (1 + X)
        kK = 4e25 * Z * (1 + X) * rho * T**(-7/2)
        kH = 1.1e-25 * np.sqrt(Z * rho) * T**7.7
        km = 0.1 * Z
        kappa = km + 1/(kH**-1 + (ke + kK)**-1)
        if log:
            return np.log10(kappa)
        else:
            return kappa
    else:
        raise Exception("Please choose a valid analytic opactiy prescription.")


# kappa_data = np.genfromtxt('X07Y029Z001.txt')
# kappa_data = np.genfromtxt('X07Y027Z003.txt')
# kappa_data = np.genfromtxt('A96X070Y027Z003.txt')
kappa_data = np.genfromtxt('TabulatedOpacities/X070Y028Z002.txt')
kappa_R = kappa_data[0, 1:]
kappa_T = kappa_data[1:, 0]
kappa_interp = interp.RegularGridInterpolator((kappa_T, kappa_R), kappa_data[1:, 1:], bounds_error=False)
kappa_data = kappa_data[1:, 1:]

lowtemp_kappa_data = np.genfromtxt('TabulatedOpacities/lowtempX07Y028Z002.txt')[::-1, :]
low_kappa_R = lowtemp_kappa_data[-1, 1:]
low_kappa_T = lowtemp_kappa_data[:-1, 0]
low_kappa_interp = interp.RegularGridInterpolator((low_kappa_T, low_kappa_R), lowtemp_kappa_data[1:, :-1], bounds_error=False)
lowtemp_kappa_data = lowtemp_kappa_data[:-1, 1:]

kappa_T_overlap = np.intersect1d(low_kappa_T, kappa_T)


def kappa_from_data(logrho, logT):
    '''
    '''
    if logT < 3:
        return 0.76
    logR = logrho - (3. * (logT - 6.))
    if logT <= low_kappa_T[-1]:
        if logR < low_kappa_R[0]:
            if logT < low_kappa_T[0]:
                return lowtemp_kappa_data[0, 0]
            else:
                return low_kappa_interp([logT, low_kappa_R[0]])
        elif logR > low_kappa_R[-1]:
            if logT < low_kappa_T[0]:
                return lowtemp_kappa_data[0, -1]
            else:
                return low_kappa_interp([logT, low_kappa_R[-1]])
        elif logT < low_kappa_T[0]:
            return low_kappa_interp([low_kappa_T[0], logR])
        elif logT < kappa_T[0] and logT > low_kappa_T[0]:
            return low_kappa_interp([logT, logR])
        prop = np.interp(logT, kappa_T_overlap, np.linspace(0, 1, len(kappa_T_overlap)))
        return prop * kappa_interp([logT, logR]) + (1 - prop) * low_kappa_interp([logT, logR])
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


def log_system(x, r, Mdot, angvel, alpha, c_cgs, m_cgs, min_Q, b):
    '''
    '''
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
    if not np.isscalar(eq11):
        eq11 = eq11[0]
    # eq11 = x[3] - kappa_from_formula(10.**x[6], 10.**x[1], log=True)
    eq12 = x[6] - (np.log10(angvel**2 / (2 * np.pi * G_cgs)) - max(x[8], np.log10(min_Q)))
    # if Q <= min_Q, density is analytic and effective temperature isn't
    if x[8] <= np.log10(min_Q):
        eq1 = 4 * x[1] - (4 * x[0] + np.log10(3 * 10.**x[2] / 8 + 0.5 + 1 / (4 * 10.**x[2])))
    return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12]


def angvel(r, M):
    '''
    '''
    return np.sqrt(G_cgs * M * M_odot_cgs / r**3)


def disk_model(M, f_edd, alpha, b):
    '''
    '''
    # Teff = x[0], T = x[1], tau = x[2], kappa = x[3], Sigma = x[4], cs = x[5], rho = x[6], h = x[7], Q = x[8], beta = x[9], prad = x[10], pgas = x[11]
    root = np.array([6, 6., 5., -1., 6., 8., -10., -2., 4., -3., 10., 4.])

    n = 100
    min_Q = 1.

    Mdot_edd = 4 * np.pi * G_cgs * M * M_odot_cgs * m_H / (0.1 * thomson_cgs * c_cgs)
    Mdot = f_edd * Mdot_edd
    rs_cm = 2 * G_cgs * M * M_odot_cgs / c_cgs**2
    r_min = 3 * rs_cm

    log_radii = np.logspace(1.1 * np.log10(r_min / rs_cm), 6, n)
    radii = log_radii * rs_cm

    temps, rho, Sigma, h, kappa, tau, Q = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    t_eff, cs, beta, prad, pgas = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)

    for i, r in enumerate(radii):
        Mdotdash = Mdot * (1 - np.sqrt(r_min / r))
        angvel_r = angvel(r, M)
        root = opt.fsolve(log_system, root, args=(r, Mdotdash, angvel_r, alpha, c_cgs, m_cgs, min_Q, b))
        t_eff[i], temps[i], tau[i], kappa[i], Sigma[i], cs[i], rho[i], h[i], Q[i], beta[i], prad[i], pgas[i] = root
        Q[i] = max(Q[i], min_Q)

    return [log_radii, t_eff, temps, tau, kappa, Sigma, cs, rho, h, Q, beta, prad, pgas]


def save_disk_model(disk_params, location='', name='', save_all=False):
    '''
    '''
    path = os.path.dirname(os.path.abspath(__file__)) + location
    if not os.path.isdir(path):
        os.mkdir(path)
    param_names = ['log_radii', 't_eff', 'temps', 'tau', 'kappa', 'Sigma', 'cs', 'rho', 'h', 'Q', 'beta', 'prad', 'pgas', 'pressure']
    filenames = [path + param + '_' + name + '.csv' for param in param_names]
    if save_all:
        indices = np.arange(0, len(param_names))
    else:   # only saves parameter arrays for those relevant to migration in the nbody code
        indices = [0, 2, 4, 5, 8, 13]
    for index in indices:
        if index == 0:
            np.savetxt(filenames[index], np.log10(disk_params[index]), delimiter=',')
        elif index == 13:
            np.savetxt(filenames[index], np.log10(10**disk_params[11] + 10**disk_params[12]), delimiter=',')
        else:
            np.savetxt(filenames[index], disk_params[index], delimiter=',')


def plot_disk_model(disk_params, axes=[], save=False, location=''):
    '''
    '''
    log_radii, t_eff, temps, tau, kappa, Sigma, cs, rho, h, Q, beta, prad, pgas = disk_params

    if len(axes) == 0:
        fig, axes = plt.subplots(nrows=6, sharex=True, figsize=(5, 10), gridspec_kw={'hspace': 0})

    axes[0].plot(log_radii, 10**temps)
    axes[1].plot(log_radii, 10**Sigma)
    axes[2].plot(log_radii, 10**h)
    axes[3].plot(log_radii, 10**kappa)
    axes[4].plot(log_radii, 10**tau)
    axes[5].plot(log_radii, 10**Q)

    for i, ax in enumerate(axes):
        ax.set(xscale='log', yscale='log')
        if i != 0:
            ax.xaxis.set_tick_params(which='both', reset=True)
    axes[0].set(ylabel='$T_{\mathrm{mid}}$ (K)')
    axes[1].set(ylabel='$\Sigma$ (g/cm$^2$)')
    axes[2].set(ylabel='$h$ ($H$/$r$)')
    axes[3].set(ylabel='$\kappa$ (cm$^2$/g)')
    axes[4].set(ylabel=r'$\tau$')
    axes[5].set(ylabel='Toomre, $Q$', xlabel='$r$/$R_s$')

    if save:
        path = os.path.dirname(os.path.abspath(__file__)) + location
        if not os.path.isdir(path):
            os.mkdir(path)
        fig.savefig(path + 'disk_model.png', dpi=400, bbox_inches='tight')


def plot_many_models():
    ''' Used to plot a range of disk models for display in the paper. Saves the images to the "Images" folder.
    '''
    # set LaTeX font for our figures
    plt.rcParams.update({"text.usetex": True})
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'cm'
    fig, axes = plt.subplots(nrows=6, sharex=True, figsize=(5, 10), gridspec_kw={'hspace': 0})
    masses = [1e6, 1e7, 1e8]
    fracs = [0.1, 0.5, 1]
    # alphas = [0.01, 0.1]
    alphas = [0.01]

    colours = ['tab:orange', 'tab:red', 'tab:purple']
    ls = ['-', '--', ':']
    lw = [1, 0.5]
    for i, M in enumerate(masses):
        for j, f_edd in enumerate(fracs):
            for k, alpha in enumerate(alphas):
                log_radii, t_eff, temps, tau, kappa, Sigma, cs, rho, h, Q, beta, prad, pgas = disk_model(M, f_edd, alpha, 0)

                axes[0].plot(log_radii, 10**temps, color=colours[i], ls=ls[j], lw=lw[k], rasterized=True)
                axes[1].plot(log_radii, 10**Sigma, color=colours[i], ls=ls[j], lw=lw[k], rasterized=True)
                axes[2].plot(log_radii, 10**h, color=colours[i], ls=ls[j], lw=lw[k], rasterized=True)
                axes[3].plot(log_radii, 10**kappa, color=colours[i], ls=ls[j], lw=lw[k], rasterized=True)
                axes[4].plot(log_radii, 10**tau, color=colours[i], ls=ls[j], lw=lw[k], rasterized=True)
                axes[5].plot(log_radii, 10**Q, color=colours[i], ls=ls[j], lw=lw[k], rasterized=True)

    axes[0].set(ylabel='$T_{\mathrm{mid}}$ (K)')
    axes[1].set(ylabel='$\Sigma$ (g/cm$^2$)')
    axes[2].set(ylabel='$h$ ($H$/$r$)')
    axes[3].set(ylabel='$\kappa$ (cm$^2$/g)')
    axes[4].set(ylabel=r'$\tau$')
    axes[5].set(ylabel='Toomre, $Q$', xlabel='$R/R_s$')
    for i, ax in enumerate(axes):
        ax.set(xscale='log', yscale='log')

    from matplotlib.lines import Line2D
    custom_lines1 = [Line2D([0], [0], color=colours[0]),
                     Line2D([0], [0], color=colours[1]),
                     Line2D([0], [0], color=colours[2])]
    custom_lines2 = [Line2D([0], [0], color='k', ls=ls[0]),
                     Line2D([0], [0], color='k', ls=ls[1]),
                     Line2D([0], [0], color='k', ls=ls[2])]
    axes[0].legend(custom_lines1, ['$M=10^6 M_\odot$',
                   '$M=10^7 M_\odot$', '$M=10^8 M_\odot$'])
    axes[-1].legend(custom_lines2, ['$f_{\mathrm{edd}} = 0.1$',
                    '$f_{\mathrm{edd}} = 0.5$', '$f_{\mathrm{edd}} = 1$'])

    fig.savefig("Images/SGDiskModels.png", dpi=400, bbox_inches='tight')
    fig.savefig("Images/SGDiskModels.pdf", dpi=400, bbox_inches='tight')

def plot_torques(M, f_edd, visc, disk_params, save=False, location=''):
    '''
    Open problems:
        1. Not sure whether to model based on total pressure or just gas pressure. Evgeni modelled by gas pressure, and this
            means that there are some migration traps in the inner disk; these migration traps disappear when modelling via total pressure
        2. Need to plot regions of parameter space that contain at least one migration trap (and at what radius!)
    '''
    import scipy.interpolate as interp

    fig, ax = plt.subplots()

    accretion = 0.1
    
    bh_mass = 10 * M_odot_cgs
    gamma_coeff = 5/3 

    rs = 2 * G_cgs * M * M_odot_cgs / c_cgs**2
    log_radii, t_eff, temps, tau, kappa, Sigma, cs, rho, h, Q, beta, prad, pgas = disk_model(M, f_edd, visc, 0)

    spl_sigma = interp.CubicSpline(np.log10(log_radii), Sigma, extrapolate=True)
    spl_temp = interp.CubicSpline(np.log10(log_radii), temps, extrapolate=True)
    spl_dens = interp.CubicSpline(np.log10(log_radii), rho, extrapolate=True)
    spl_h = interp.CubicSpline(np.log10(log_radii), h, extrapolate=True)
    spl_kappa = interp.CubicSpline(np.log10(log_radii), kappa, extrapolate=True)
    spl_tau = interp.CubicSpline(np.log10(log_radii), tau, extrapolate=True)
    spl_P = interp.CubicSpline(np.log10(log_radii), np.log10(10**prad + 10**pgas), extrapolate=True) # total pressure
    # spl_P = interp.CubicSpline(np.log10(log_radii), pgas, extrapolate=True)     # radiation pressure
    spl_cs = interp.CubicSpline(np.log10(log_radii), cs, extrapolate=True)

    def alpha(r): return -spl_sigma.derivative()(np.log10(r))
    def beta(r): return -spl_temp.derivative()(np.log10(r))
    def P_deriv(r): return -spl_P.derivative()(np.log10(r))
    
    log_radii = np.logspace(1, 5, 1000)
    torques = np.zeros(len(log_radii))

    for ii, r in enumerate(log_radii):
        logr = np.log10(r)
        Gamma_0 = ((10/M) / 10**spl_h(logr))**2 * 10**spl_sigma(logr) * (r*rs)**4 * angvel(r*rs, M)**2
        
        ### Migration from pardekooper
        # c_v = 14304 / 1000
        # tau_eff = 3 * 10**spl_tau(logr) / 8 + np.sqrt(3)/4 + 0.25 / 10**spl_tau(logr)
        # Theta = (c_v * 10**spl_sigma(logr) * angvel(r*rs, M) * tau_eff) / (12. * np.pi * stef_boltz * 10**(3 * spl_temp(logr)));
        # Gamma_iso = -0.85 - alpha(r) - 0.9 * beta(r)
        # xi = beta(r) - (gamma_coeff - 1) * alpha(r)
        # Gamma_ad = (-0.85 - alpha(r) - 1.7 * beta(r) + 7.9 * xi / gamma_coeff) / gamma_coeff;
        # Gamma = Gamma_0 * (Gamma_ad * Theta*Theta + Gamma_iso) / ((Theta + 1)*(Theta + 1));
        
        ### Migration from Jimenez
        H = 10**spl_h(logr) * r*rs
        chi = 16. * gamma_coeff * (gamma_coeff - 1.) * stef_boltz * 10**(4 * spl_temp(logr)) / (3. * 10**(2 * spl_dens(logr)) * 10**spl_kappa(logr) * (angvel(r*rs, M) * H)**2)
        chi_chi_c = chi / (H**2 * angvel(r*rs, M))
        fx = (np.sqrt(chi_chi_c / 2.) + 1. / gamma_coeff) / (np.sqrt(chi_chi_c / 2.) + 1.);
        Gamma_lindblad = - (2.34 - 0.1 * alpha(r) + 1.5 * beta(r)) * fx;
        Gamma_simp_corot = (0.46 - 0.96 * alpha(r) + 1.8 * beta(r)) / gamma_coeff;
        Gamma = Gamma_0 * (Gamma_lindblad + Gamma_simp_corot);

        ### Thermal torques
        dPdr = P_deriv(r)
        x_c = dPdr * H**2 / (3 * gamma_coeff * r*rs)
        L = accretion * 4. * np.pi * G_cgs * bh_mass * m_H * c_cgs / thomson_cgs;     # accretion assuming completely ionized hydrogen
        # L = accretion * 4 * np.pi * G_cgs * bh_mass * c_cgs / 10**spl_kappa(logr)       # accretion assuming the AGN disk composition
        # below are equations 17-23 from gilbaum 2022
        R_BHL = 2 * G_cgs * bh_mass / (H * angvel(r*rs, M))**2
        R_H = r*rs * np.cbrt(10 / (3 * M))
        b_H = np.sqrt(R_BHL * R_H)
        mdot_RBHL = np.pi * min(R_BHL, b_H) * min(R_BHL, b_H, H) * (H * angvel(r*rs, M))
        L_RBHL = 0.1 * c_cgs**2 * mdot_RBHL
        L = min(L_RBHL, L / accretion)
        
        Lc = 4. * np.pi * G_cgs * bh_mass * 10**spl_dens(logr) * chi / gamma_coeff
        # print(L/Lc)
        lambda_ = np.sqrt(2. * chi / (3 * gamma_coeff * angvel(r*rs, M)));
        Gamma_thermal = 1.61 * (gamma_coeff - 1) / gamma_coeff * x_c / lambda_ * (L/Lc - 1.) * Gamma_0 / 10**spl_h(logr);

        ### GR Inspiral torque
        Gamma_GW = Gamma_0 * (-32 / 5 * (c_cgs / 10**spl_cs(logr))**3 * 10**(6 * spl_h(logr)) * (2*r)**-4 * M*M_odot_cgs / (10**spl_sigma(logr) * (r*rs)**2))
        
        Gamma += Gamma_thermal + Gamma_GW
        torques[ii] = Gamma
        
    pos_vals = torques > 0 
    neg_vals = torques <= 0 
    pos_torques = [torques[i] if pos_vals[i] else np.nan for i in range(len(torques))]
    neg_torques = [-torques[i] if neg_vals[i] else np.nan for i in range(len(torques))]
    ax.plot(log_radii, pos_torques, label='$+$ve', rasterized=True)
    ax.plot(log_radii, neg_torques, ls='--', label='$-$ve', rasterized=True)
    
    ax.set(xscale='log', yscale='log', xlabel="log$(R/R_s)$", ylabel='abs($\Gamma$)')
    ax.legend()
    ax.grid()
    if save:
        path = os.path.dirname(os.path.abspath(__file__)) + location
        if not os.path.isdir(path):
            os.mkdir(path)
        fig.savefig(path + 'Torque_Model.png', dpi=400, bbox_inches='tight')
        fig.savefig(path + 'Torque_Model.pdf', dpi=400, bbox_inches='tight')

def plot_many_torques():
    '''
    Open problems:
        1. Not sure whether to model based on total pressure or just gas pressure. Evgeni modelled by gas pressure, and this
            means that there are some migration traps in the inner disk; these migration traps disappear when modelling via total pressure
        2. Need to plot regions of parameter space that contain at least one migration trap (and at what radius!)
    '''
    import scipy.interpolate as interp
    plt.rcParams.update({"text.usetex": True})
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'cm'

    fig, ax = plt.subplots()
    # fig2, ax2 = plt.subplots()

    accretion = 1
    masses = [1e6, 1e7, 1e8, 1e9]
    fracs = [0.1]
    # alphas = [0.01, 0.1]
    alphas = [0.01]
    colours = ['tab:orange', 'tab:red', 'tab:purple', 'tab:blue']
    ls = ['-', '--', ':']
    lw = [1, 0.5]
    
    bh_mass = 10 * M_odot_cgs
    # R_mu = 8.3145 / 2.016
    R_mu = 8.3145 / (2.016 * m_H * 6.022e23 * 10)
    gamma_coeff = 5/3 

    for i, M in enumerate(masses):
        for j, f_edd in enumerate(fracs):
            for jj, visc in enumerate(alphas):
                rs = 2 * G_cgs * M * M_odot_cgs / c_cgs**2
                log_radii, t_eff, temps, tau, kappa, Sigma, cs, rho, h, Q, beta, prad, pgas = disk_model(M, f_edd, visc, 0)

                spl_sigma = interp.CubicSpline(np.log10(log_radii), Sigma, extrapolate=True)
                spl_temp = interp.CubicSpline(np.log10(log_radii), temps, extrapolate=True)
                spl_dens = interp.CubicSpline(np.log10(log_radii), rho, extrapolate=True)
                spl_h = interp.CubicSpline(np.log10(log_radii), h, extrapolate=True)
                spl_kappa = interp.CubicSpline(np.log10(log_radii), kappa, extrapolate=True)
                spl_tau = interp.CubicSpline(np.log10(log_radii), tau, extrapolate=True)
                spl_P = interp.CubicSpline(np.log10(log_radii), np.log10(10**prad + 10**pgas), extrapolate=True) # total pressure
                # spl_P = interp.CubicSpline(np.log10(log_radii), pgas, extrapolate=True)     # radiation pressure
                spl_cs = interp.CubicSpline(np.log10(log_radii), cs, extrapolate=True)

                def alpha(r): return -spl_sigma.derivative()(np.log10(r))
                def beta(r): return -spl_temp.derivative()(np.log10(r))
                def P_deriv(r): return -spl_P.derivative()(np.log10(r))
                
                log_radii = np.logspace(1, 5, 1000)
                torques = np.zeros(len(log_radii))

                for ii, r in enumerate(log_radii):
                    logr = np.log10(r)
                    Gamma_0 = ((10/M) / 10**spl_h(logr))**2 * 10**spl_sigma(logr) * (r*rs)**4 * angvel(r*rs, M)**2
                    
                    ### Migration from pardekooper
                    # c_v = 14304 / 1000
                    # tau_eff = 3 * 10**spl_tau(logr) / 8 + np.sqrt(3)/4 + 0.25 / 10**spl_tau(logr)
                    # Theta = (c_v * 10**spl_sigma(logr) * angvel(r*rs, M) * tau_eff) / (12. * np.pi * stef_boltz * 10**(3 * spl_temp(logr)));
                    # Gamma_iso = -0.85 - alpha(r) - 0.9 * beta(r)
                    # xi = beta(r) - (gamma_coeff - 1) * alpha(r)
                    # Gamma_ad = (-0.85 - alpha(r) - 1.7 * beta(r) + 7.9 * xi / gamma_coeff) / gamma_coeff;
                    # Gamma = Gamma_0 * (Gamma_ad * Theta*Theta + Gamma_iso) / ((Theta + 1)*(Theta + 1));
                    
                    ### Migration from Jimenez
                    cs = 10**spl_cs(logr)
                    H = 10**spl_h(logr) * r*rs
                    chi = 16. * gamma_coeff * (gamma_coeff - 1.) * stef_boltz * 10**(4 * spl_temp(logr)) / (3. * 10**(2 * spl_dens(logr)) * 10**spl_kappa(logr) * cs**2)
                    chi_chi_c = chi / (H**2 * angvel(r*rs, M))
                    fx = (np.sqrt(chi_chi_c / 2.) + 1. / gamma_coeff) / (np.sqrt(chi_chi_c / 2.) + 1.);
                    Gamma_lindblad = - (2.34 - 0.1 * alpha(r) + 1.5 * beta(r)) * fx;
                    Gamma_simp_corot = (0.46 - 0.96 * alpha(r) + 1.8 * beta(r)) / gamma_coeff;
                    Gamma = Gamma_0 * (Gamma_lindblad + Gamma_simp_corot);
                    
                    ## Thermal torques
                    dPdr = P_deriv(r)
                    x_c = dPdr * H**2 / (3 * gamma_coeff * r*rs)
                    L = accretion * 4. * np.pi * G_cgs * bh_mass * m_H * c_cgs / thomson_cgs;     # accretion assuming completely ionized hydrogen
                    # below are equations 17-23 from gilbaum 2022
                    R_BHL = 2 * G_cgs * bh_mass / cs**2
                    R_H = r*rs * np.cbrt(10 / (3 * M))
                    b_H = np.sqrt(R_BHL * R_H)
                    mdot_RBHL = np.pi * min(R_BHL, b_H) * min(R_BHL, b_H, H) * cs
                    L_RBHL = 0.1 * c_cgs**2 * mdot_RBHL
                    L = min(L_RBHL, L / accretion)
                    Lc = 4. * np.pi * G_cgs * bh_mass * 10**spl_dens(logr) * chi / gamma_coeff
                    lambda_ = np.sqrt(2. * chi / (3 * gamma_coeff * angvel(r*rs, M)));
                    Gamma_thermal = 1.61 * (gamma_coeff - 1) / gamma_coeff * x_c / lambda_ * (L/Lc - 1.) * Gamma_0 / 10**spl_h(logr);
                    
                    ### GR Inspiral torque
                    Gamma_GW = Gamma_0 * (-32 / 5 * (c_cgs / cs)**3 * 10**(6 * spl_h(logr)) * (2*r)**-4 * M*M_odot_cgs / (10**spl_sigma(logr) * (r*rs)**2))
                    # if M==1e8 and (1e3 <= r <= 4e3): 
                        # print(x_c / lambda_)
                        # print(Gamma_thermal / Gamma)
                        # print(L/Lc)
                        # print(chi_chi_c)
                    #     ax2.scatter(r, Gamma_thermal/Gamma, s=1)
                    Gamma += Gamma_thermal + Gamma_GW
                    torques[ii] = Gamma
                    
                    
                pos_vals = torques > 0 
                neg_vals = torques <= 0 
                pos_torques = [torques[i] if pos_vals[i] else np.nan for i in range(len(torques))]
                neg_torques = [-torques[i] if neg_vals[i] else np.nan for i in range(len(torques))]
                ax.plot(log_radii, pos_torques, c=colours[i], label=f'$M=10^{int(np.log10(M))}$', rasterized=True)
                ax.plot(log_radii, neg_torques, c=colours[i], ls='--', rasterized=True)
    ax.set(xscale='log', yscale='log', xlabel="log$(R/R_s)$", ylabel='abs($\Gamma$)')
    ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    from matplotlib.lines import Line2D
    p1 = Line2D([0], [0], color='k', ls='-', label='$+$ve'); handles.append(p1)
    p2 = Line2D([0], [0], color='k', ls='--', label='$-$ve'); handles.append(p2)
    ax.legend(handles=handles)
    ax.grid()
    fig.savefig('Images/Torque_Model.png', dpi=400, bbox_inches='tight')
    fig.savefig('Images/Torque_Model.pdf', dpi=400, bbox_inches='tight')
    
    # ax2.set(xscale='log')

def plot_migration_traps(visc):
    fig, ax = plt.subplots()
    n_M = 20
    n_edd = 20
    bh_mass = 10 * M_odot_cgs
    gamma_coeff = 5/3 
    trap_rads = np.ones((n_M, n_edd))
    masses = np.logspace(6, 9, n_M)
    fractions = np.logspace(-2, 0, n_edd)
    for i, M in enumerate(masses):
        for j, f_edd in enumerate(fractions):
            rs = 2 * G_cgs * M * M_odot_cgs / c_cgs**2
            log_radii, t_eff, temps, tau, kappa, Sigma, cs, rho, h, Q, beta, prad, pgas = disk_model(M, f_edd, visc, 0)
    
            spl_sigma = interp.CubicSpline(np.log10(log_radii), Sigma, extrapolate=True)
            spl_temp = interp.CubicSpline(np.log10(log_radii), temps, extrapolate=True)
            spl_dens = interp.CubicSpline(np.log10(log_radii), rho, extrapolate=True)
            spl_h = interp.CubicSpline(np.log10(log_radii), h, extrapolate=True)
            spl_kappa = interp.CubicSpline(np.log10(log_radii), kappa, extrapolate=True)
            spl_tau = interp.CubicSpline(np.log10(log_radii), tau, extrapolate=True)
            spl_P = interp.CubicSpline(np.log10(log_radii), np.log10(10**prad + 10**pgas), extrapolate=True) # total pressure
            spl_cs = interp.CubicSpline(np.log10(log_radii), cs, extrapolate=True)
    
            def alpha(r): return -spl_sigma.derivative()(np.log10(r))
            def beta(r): return -spl_temp.derivative()(np.log10(r))
            def P_deriv(r): return -spl_P.derivative()(np.log10(r))
            
            log_radii = np.logspace(1, 5, 1000)
            torques = np.zeros(len(log_radii))
    
            for ii, r in enumerate(log_radii):
                logr = np.log10(r)
                Gamma_0 = ((10/M) / 10**spl_h(logr))**2 * 10**spl_sigma(logr) * (r*rs)**4 * angvel(r*rs, M)**2
                
                ### Migration from pardekooper
                # c_v = 14304 / 1000
                # tau_eff = 3 * 10**spl_tau(logr) / 8 + np.sqrt(3)/4 + 0.25 / 10**spl_tau(logr)
                # Theta = (c_v * 10**spl_sigma(logr) * angvel(r*rs, M) * tau_eff) / (12. * np.pi * stef_boltz * 10**(3 * spl_temp(logr)));
                # Gamma_iso = -0.85 - alpha(r) - 0.9 * beta(r)
                # xi = beta(r) - (gamma_coeff - 1) * alpha(r)
                # Gamma_ad = (-0.85 - alpha(r) - 1.7 * beta(r) + 7.9 * xi / gamma_coeff) / gamma_coeff;
                # Gamma = Gamma_0 * (Gamma_ad * Theta*Theta + Gamma_iso) / ((Theta + 1)*(Theta + 1));
                
                ### Migration from Jimenez
                cs = 10**spl_cs(logr)
                H = 10**spl_h(logr) * r*rs
                chi = 16. * gamma_coeff * (gamma_coeff - 1.) * stef_boltz * 10**(4 * spl_temp(logr)) / (3. * 10**(2 * spl_dens(logr)) * 10**spl_kappa(logr) * cs**2)
                chi_chi_c = chi / (H**2 * angvel(r*rs, M))
                fx = (np.sqrt(chi_chi_c / 2.) + 1. / gamma_coeff) / (np.sqrt(chi_chi_c / 2.) + 1.);
                Gamma_lindblad = - (2.34 - 0.1 * alpha(r) + 1.5 * beta(r)) * fx;
                Gamma_simp_corot = (0.46 - 0.96 * alpha(r) + 1.8 * beta(r)) / gamma_coeff;
                Gamma = Gamma_0 * (Gamma_lindblad + Gamma_simp_corot);
                
                ## Thermal torques
                dPdr = P_deriv(r)
                x_c = dPdr * H**2 / (3 * gamma_coeff * r*rs)
                L = 4. * np.pi * G_cgs * bh_mass * m_H * c_cgs / thomson_cgs;     # accretion assuming completely ionized hydrogen
                # below are equations 17-23 from gilbaum 2022
                # R_BHL = 2 * G_cgs * bh_mass / cs**2
                # R_H = r*rs * np.cbrt(10 / (3 * M))
                # b_H = np.sqrt(R_BHL * R_H)
                # mdot_RBHL = np.pi * min(R_BHL, b_H) * min(R_BHL, b_H, H) * cs
                # L_RBHL = 0.1 * c_cgs**2 * mdot_RBHL
                # L = min(L_RBHL, L)
                Lc = 4. * np.pi * G_cgs * bh_mass * 10**spl_dens(logr) * chi / gamma_coeff
                lambda_ = np.sqrt(2. * chi / (3 * gamma_coeff * angvel(r*rs, M)));
                Gamma_thermal = 1.61 * (gamma_coeff - 1) / gamma_coeff * x_c / lambda_ * (L/Lc - 1.) * Gamma_0 / 10**spl_h(logr);
                
                ### GR Inspiral torque
                Gamma_GW = Gamma_0 * (-32 / 5 * (c_cgs / cs)**3 * 10**(6 * spl_h(logr)) * (2*r)**-4 * M*M_odot_cgs / (10**spl_sigma(logr) * (r*rs)**2))
    
                Gamma += Gamma_thermal + Gamma_GW
                torques[ii] = Gamma
                
            if torques[-1] < 0:
                for n_torque, torque in enumerate(torques[::-1]):
                    k = len(torques) - n_torque
                    if torque >= 0:
                        trap_rads[i, j] = log_radii[k]
                        break
            # fig2, ax2 = plt.subplots()
            # pos_vals = torques > 0 
            # neg_vals = torques <= 0 
            # pos_torques = [torques[iii] if pos_vals[iii] else np.nan for iii in range(len(torques))]
            # neg_torques = [-torques[iii] if neg_vals[iii] else np.nan for iii in range(len(torques))]
            # ax2.plot(log_radii, pos_torques, rasterized=True)
            # ax2.plot(log_radii, neg_torques, ls='--', rasterized=True)
            # ax2.set(xscale='log', yscale='log')
    x, y = np.meshgrid(masses, fractions)
    from matplotlib.colors import LogNorm
    # from matplotlib import cm, ticker
    contour = ax.pcolormesh(x, y, trap_rads, cmap='viridis', 
                            norm=LogNorm(vmin=trap_rads.min(), vmax=trap_rads.max()),
                            rasterized=True)
    ax.set(xlabel='SMBH Mass', ylabel='Eddington Fraction', xscale='log', yscale='log')
    fig.colorbar(contour, label='Migration Trap Location ($R_s$)')
    fig.savefig(f'Images/MigrationTraps-alph{visc}.png', dpi=400, bbox_inches='tight')
    fig.savefig(f'Images/MigrationTraps-alph{visc}.pdf', dpi=400, bbox_inches='tight')
                
                    
            

# plot_many_models()
# plot_many_torques()
# plot_migration_traps(0.01)