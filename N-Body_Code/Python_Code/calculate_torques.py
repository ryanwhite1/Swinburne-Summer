#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 16:39:12 2023

@author: evgeni
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# plt.rc('text', usetex=False)
# matplotlib.rcParams.update({'font.size': 20})

#constants in cgs
G=6.674e-8
msun = 1.989e33
au = 1.49e13
c=2.99792458e10
sb=5.6e-5

#log space of spatial separations, in units of r_g = GM_smbh/c^2
rs= np.logspace(0.01,6,1000)

### construct a locat solution

### Input
# r = location in units of rg
# m = mass of SMBH in units of 1e8 M_\odot
# m_dot = accretion rate in units of the critical accretion rate (see Shakura & Sunyaev 1973)
# alpha = the alpha viscosity
# m_stellar= the mass of the stella BH, in units of M_\odot

### returns a vector
#  rho, H, cs, P, Sigma, T, kappa, zone, kappa_m17, P_grad, Sigma_grad, T_grad, gamma, Gamma_0
# rho = density
# H = scale height 
# P pressure i nthe midplane
# Sigma = surface density
# T = temperature at the midplane
# kappa = opacity as power laws in dfferent zones
# zone = the zone for which the solution is valid
#kappa_m17 = the opacity from Metzger and Pejcha (2017)
# P_ grad = - d lnP / d lnr = pressure density gradient
# Sigma_grad = - d lnSigma / d lnr = surface density gradient
# T_grad = - d lnT/ d lnr = temperature gradient
# gamma = adiabatic index 
# Gamma_0 = the nominal torque at each location
def get_disc_params(r,m, m_dot, alpha, m_stellar):
    # transition radii to the different zones
    R12 = 449.842 * alpha**(2/21) * m**(2/21) * m_dot**(16/21)
    R23 = 987.891 * m_dot**(2/3)
#    R24 =  1383.38 * alpha**-0.24 * m**-0.24 * m_dot**0.423
    R34 = 3333.58 * alpha**(-0.28) * m**(-0.28) * m_dot**0.385
    R1Q = 498 * alpha**(2/9) * m**(-2/9) * m_dot **(4/9)
    R2Q = 634.4 * alpha**(14/27) * m**(-26/27) * m_dot**(-8/27)
    R3Q = 580.65 * alpha**(28/45) * m**(-52/45) * m_dot**(-22/45)
    R4Q = 184.709 * alpha**(14/27) * m**(-26/27) * m_dot**(-8/27)

    X=0.7381
    Y=0.2485
    Z=0.0134
    kappa_es=0.2*(1+X)
    m_dot = m_dot * (1-r**-0.5)
    
    # zone 1: P = P_rad; kappa = kappa_es (Shakura & Sunyaev 1973)
    if r<=R12 and r<=R1Q:
        zone=1
        gamma=4/3
        rho = 3.29e-14* alpha **-1 * m**-1 * m_dot**-2 * r**1.5 
        H = 1.3e14 * m * m_dot
        cs = 1.8e10 * m_dot * r**-1.5
        P = 1.065e7 * alpha**-1 * m**-1 * r**-1.5
        Sigma = 8.574 * alpha**-1 * m_dot**-1 * r**1.5 
        T = 254920 * alpha**-0.25 * m**-0.25 * r**(-3/8)
        kappa=kappa_es
    
        kappa_R = 4e25*(1+X)*Z * rho * T**-3.5
        kappa_H_minus = 1.1e-25*Z**0.5 * rho**0.5 * T**7.7
        kappa_m17 = (kappa_H_minus**-1 + (kappa_es + kappa_R)**-1)**-1
        
        P_grad=1.5
        Sigma_grad=-1.5
        T_grad = 3/8

    # zone 2: P = P_gas; kappa = kappa_es (Shakura & Sunyaev 1973)
    # the transition to gas pressure usually occurs in a hot environment wher the opacirt still due to ES
    elif R12 <= r <= R23 and r<=R2Q:
        zone=2
        gamma=5/3
        rho = 7.495e-6 * alpha**-0.7 * m**-0.7 * m_dot**(2/5) * r**(-33/20)
        H = 2.13e11 * m**0.9 * m_dot**0.2 * r**(21/20) * alpha**-0.1
        cs = 2.945e7 * m**-0.1 *alpha**-0.1 * r**(-9/20) * m_dot**0.2
        P = 6.5e9 * alpha**-0.9 * m**-0.9 * r**(-51/20) * m_dot**0.8
        Sigma = 3.196e6 * alpha**-0.8 * m**0.2 * m_dot**(3/5) * r**(-3/5)
        T = 6.3e6 * alpha**-0.2 * m**-0.2 * m_dot**(2/5) * r**(-9/10)
        kappa=kappa_es
        
        kappa_R = 4e25*(1+X)*Z * rho * T**-3.5
        kappa_H_minus = 1.1e-25*Z**0.5 * rho**0.5 * T**7.7
        kappa_m17 = (kappa_H_minus**-1 + (kappa_es + kappa_R)**-1)**-1

        P_grad=51/20
        Sigma_grad=3/5
        T_grad = 9/10

    # zone 3: P = P_gas; kappa = kappa_kramers (Shakura & Sunyaev 1973)
    # 
    elif R23 <=r <= R34 and r<=R3Q:
        zone=3
        gamma=5/3
        rho =  3.536e-5* alpha**-0.7 * m**-0.7 * m_dot**(11/20) * r**(-15/8)
        H = 1.27e11 * alpha**-0.1 * m**0.9 * m_dot**(3/20) * r**(9/8)
        cs = 1.756e7 * alpha**-0.1 * m**-0.1 * m_dot**(3/20) * r**(-3/8)
        P = 1.09e10 * alpha**-0.9 * m**-0.9 * m_dot**(17/20) * r**(-21/8)
        Sigma = 9e6 * alpha**-0.8 * m**0.2 * m_dot**0.7 * r**(-3/4)
        T = 2.239e6 * alpha**-0.2* m**-0.2 * m_dot**0.3 * r**(-3/4)
        kappa_R = 4e25*(1+X)*Z * rho * T**-3.5
        kappa_H_minus = 1.1e-25*Z**0.5 * rho**0.5 * T**7.7
        
        kappa_m17 = (kappa_H_minus**-1 + (kappa_es + kappa_R)**-1)**-1
        kappa = kappa_R

        P_grad=21/8
        Sigma_grad=3/4
        T_grad = 3/4
        
        # zone 4: P = P_gas; kappa = molecular 
        # We tried to implement a narrow H-minus zone but it didn't work out. 
        # the condition to transition to this zone (i.e. R34) are indirectly determined via kappa_H_minus = kappa_molecular
        # this causes a jump in the solution at this transition, since the opacity jump is discontinuous
        # we tested against OPAL tables and comfirmed that the jump is real and there's no H_minus zone in these solution either
    elif R34 <=r <=R4Q: 
        zone=4
        gamma=5/3
        rho = 3.96e-5 * alpha**-0.7 * m**-0.7 * m_dot**(2/5) * r**(-33/20)
        H = 1.22e11 * m**0.9 * m_dot**0.2 * r**(21/20) * alpha**-0.1
        cs = 1.69e7 * m**-0.1 *alpha**-0.1 * r**(-9/20) * m_dot**0.2
        P = 1.13e10 * alpha**-0.9 * m**-0.9 * r**(-51/20) * m_dot**0.8
        Sigma = 9.7e6 * alpha**-0.8 * m**0.2 * m_dot**(3/5) * r**(-3/5)
        T = 2.07e6 * alpha**-0.2 * m**-0.2 * m_dot**(2/5) * r**(-9/10)
        kappa= 0.1 * Z
        kappa_m17 = kappa
        
        P_grad=51/20
        Sigma_grad=3/5
        T_grad = 9/10
     
        P_grad=4.54
        Sigma_grad=4.58
        T_grad = -3.08

    #Q_T = 1 : Transition to gravitationally unstable zone where Q_toomre = 1
    else:
        zone=5
        rho = 4.54949e-2 * m **-2 * r**-3
        H = 1.169e10 * alpha**(-1/3) * m**(4/3) * m_dot**(1/3) * r**1.5
        cs = 1.6146e6 *alpha**(-1/3) * m**(1/3) * m_dot**(1/3) 
        P = 1.18e11 * alpha**(-2/3) * m**(-4/3) * m_dot**(2/3) * r**-3
        Sigma = 1.064e9 * alpha**(-1/3) * m**(-2/3) * m_dot**(1/3) * r**-1.5

        # here we need to do some tricks because we don't know if the pressure is rad or gas
        # let's define some more constants in cgs        
        mu = 0.62; # mean molecular weight
        kb=1.38e-16; # Boltzmann constant 
        mp=1.67e-24 # proton mass
        
        # a critical temperature where the gas and rad pressure are equal
        T_0 = (3*c*kb/4/sb/mu/mp*rho)**(1/3)
        # temperature from gas pressure
        T_g = P/rho*mu*mp/kb;
        # temperature from radiation pressure
        T_rad = (3*c*P/4/sb)**0.25
        # assume first we're in the gas pressure regime
        T = T_g
        gamma=5/3
        T_grad = 0
        # check if we're actually in the radiation pressure
        if T_rad> T_0:
            T = T_rad
            gamma=4/3
            T_grad = 3/4
            #furthermore, T is flat in the optically thin limit
        T = T / min(1,0.1*Z*Sigma/2)**0.5
        
        kappa_R=4e25*(1+X)*Z * rho * T**-3.5
        kappa_H_minus =  1.1e-25*Z**0.5 * rho**0.5 * T**7.7
        kappa =  max(0.1*Z*1, (kappa_H_minus**-1 + (kappa_es + kappa_R)**-1)**-1)

        kappa_m17 = max(0.1*Z, (kappa_H_minus**-1 + (kappa_es + kappa_R)**-1)**-1)
  
        P_grad=3
        Sigma_grad=1.5
        #optically thin limit
        if 0.1*Z*Sigma/2 <= 1:
            T_grad=0
    
    rg = G*msun*1e8*m/c**2 # gravitational radius in CGS        
    Gamma_0 = Sigma*(r*rg/H)**5 * (r*rg)**2 * cs**2 * (m_stellar/1e8/m)**2  
    return rho, H, cs, P, Sigma, T, kappa, zone, kappa_m17, P_grad, Sigma_grad, T_grad, gamma, Gamma_0

### calcualte additional quantities
# input:
# mm, m_stellar, m_dot, alpha same as above
# which_prefactor = a string that determines which law of type I migration to use
# ##possible values are, uncomment one
#### the formulae used in Paardekooper 2010
#which_prefactor = 'p10' 
####  The furmula used in Jimenez and Masset (2017) for the isothermal case
#which_prefactor = 'JM_iso'
#### The formula used in Jimenez and Masset (2017) for the total linear case
#which_prefactor = 'JM_tot'
    
# additional quantities that can be derived. Returns a vector
# chi, lambda, x_cs, r_Hill, Gamma_0, Gamma_I, l_ratiom Gamma_thermal
# chi = thermal conductivity
# lambda = thermal diffusion length. It should be sandwicked between x_c and H for thermal torque to be efficient.
# x_cs = the corotation radius = where the tangential velocity of gas equals the Keplerian velocity of the stellar BH
# (The locations is slightly different due to the pressure gradient lowering which affects the gas motion)
# r_Hill = the Hill radius

def get_disc_torques(mm,m_stellar, m_dot,alpha, which_prefactor):
    # gravitational radius
    rg = G*msun*1e8*mm/c/c
 #  R12 = 449.842 * alpha**(2/21) * mm**(2/21) * m_dot**(16/21)

    # get the disc structure as array for all lengths    
    rhos, Hs, css, Ps, Sigmas, Ts, kappas, zoness, kappa_m17s, P_grad, Sigma_grad, T_grad, gammas, Gamma_0 = [[get_disc_params(x,mm,m_dot,alpha, m_stellar)[i] for x in rs] for i in range(0,14)]
    taus = [x*y/2 for (x,y) in zip(kappas, Sigmas)]
#    taus_17 = [x*y/2 for (x,y) in zip(kappa_m17s, Sigmas)]

    # calculate the thermal conductivity from radiative diffusion source as array
    chis = [ 16 * gamma * (gamma-1) * sb * t**4 / 3 / k/ rho**2 /cs**2 for (t,k,cs,rho, gamma) in zip(Ts, kappas,css, rhos, gammas)]  
    # get other relevant length scales
    lambdas = [(2 * chi / 3/gamma/cs*h)**0.5 for (chi,cs,h, gamma) in zip(chis, css, Hs, gammas)]
    x_cs = [np.fabs(-P_g * h**2/gamma /r/rg/3) for (P_g,h,r, gamma) in zip(P_grad, Hs, rs, gammas)]
    r_Hills = [r*rg * (m_stellar/1e8/mm/3)**(1/3) for r in rs]

    # set the Type I migration torque depending on the prefactor law we choose
    if which_prefactor=='p10':
        Gamma_I = [(-0.8 - t_g - 0.9*sigma_g) * h/r/rg for (t_g, sigma_g, r,h) in zip(T_grad, Sigma_grad, rs, Hs)]
    if which_prefactor=='JM_iso':
        Gamma_I = [-(1.36 + 0.54*sigma_g + 0.5*t_g) * h/r/rg for (t_g, sigma_g, r,h) in zip(T_grad, Sigma_grad, rs, Hs)]
    if which_prefactor=='JM_tot':
       Gamma_I = [(- (2.34 - 0.1*sigma_g + 1.5*t_g) / gamma  + (0.46 - 0.96*sigma_g + 1.8*t_g) / gamma) * h/r/rg for (t_g, sigma_g, r,h, gamma) in zip(T_grad, Sigma_grad, rs, Hs, gammas)]      
    # find the ratio of the stellar BH luminosity (=Eddington luminosity) to the critical luminosity, L_c 
    l_ratio = [gamma*c/k/rho/chi for (r,k,rho,chi, gamma) in zip(rs,kappas,rhos, chis, gammas)]
    # calculate the total thermal torque
    Gamma_thermal = [(1-np.exp(-tt*alpha**0.5))*1.61*(gamma-1)/gamma * x/lambda_c*(l-1) for (x,lambda_c,l,tt, gamma) in zip(x_cs, lambdas, l_ratio, taus, gammas)]
    return chis, lambdas, x_cs, r_Hills, Gamma_I, l_ratio, Gamma_thermal

# disc model parameters  
# accretion rate in units of m_crit (See Shakura & sunyaev for definetion)
m_dot=0.1 
#alpha viscosity
alpha=0.1
# SMBH mass in units of 1e8 M_\odot
m_smbh = 10
# mass of stellar mass BH - in m_\odot
m_stellar = 10
# whici law for type I to use
which_prefactor = 'JM_tot' 

rhos, Hs, css, Ps, Sigmas, Ts, kappas, zoness, kappa_m17s, P_grad, Sigma_grad, T_grad, gammas, Gamma_0 = [[get_disc_params(x,m_smbh,m_dot,alpha, m_stellar)[i] for x in rs] for i in range(0,14)]
chis, lambdas, x_cs, r_Hills, Gamma_I, l_ratio, Gamma_thermal = get_disc_torques(m_smbh, m_stellar, m_dot,alpha, which_prefactor)

plt.plot()
plt.subplot(221)
plt.plot(np.log10(rs), np.log10(rhos))
plt.ylabel(r'$\rho$')
plt.subplot(222)
plt.plot(np.log10(rs),np.log10( kappas))
plt.ylabel(r'$\kappa$')
plt.subplot(223)
plt.plot(np.log10(rs), np.log10([-x for x in Gamma_thermal]))
plt.plot(np.log10(rs), np.log10([x for x in Gamma_thermal]), linestyle='dashed')
plt.ylabel(r'$\Gamma_{\rm th} / \Gamma_0$')
plt.subplot(224)
plt.plot(np.log10(rs), np.log10([-x for x in Gamma_I]))
plt.ylabel(r'$\Gamma_I / \Gamma_0$')
