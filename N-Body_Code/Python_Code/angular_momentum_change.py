# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 21:12:00 2024

@author: ryanw
"""

import numpy as np


c = 299792458
M_odot = 1.98e30
G = 6.67e-11
m1 = 50
m2 = 50
r_s = 2 * G * (m1 + m2) * M_odot / c**2
r_isco = 5 * r_s
M = 1e8
R_s = 2 * G * M * M_odot / c**2
R = 5e3 * R_s
mutual_hill = (2*R * np.cbrt((m1 + m2) / (3 * M))) / (2 - 0.333 * np.cbrt((m1 + m2) / (3 * M)))

r1f = np.array([0, 0, 1]) * r_isco
r2f = np.array([0, 0, -1]) * r_isco
r1i = r1f / r_isco * mutual_hill / 2
r2i = r2f / r_isco * mutual_hill / 2

v1f = np.sqrt(G * (m1 + m2) * M_odot / r_isco) * np.array([0, -1, 0])
v2f = np.sqrt(G * (m1 + m2) * M_odot / r_isco) * np.array([0, 1, 0])
v1i = np.sqrt(G * (m1 + m2) * M_odot / (mutual_hill / 2)) * np.array([0, -1, 0])
v2i = np.sqrt(G * (m1 + m2) * M_odot / (mutual_hill / 2)) * np.array([0, 1, 0])

Li = np.cross(r1i, m1 * M_odot * v1i) + np.cross(r2i, m2 * M_odot * v2i)
Lf = np.cross(r1f, m1 * M_odot * v1f) + np.cross(r2f, m2 * M_odot * v2f)

del_L = np.linalg.norm(Lf) / np.linalg.norm(Li)

