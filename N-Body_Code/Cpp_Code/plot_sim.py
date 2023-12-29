import numpy as np
import matplotlib.pyplot as plt

def eccentricity(position, velocity):
    '''https://en.wikipedia.org/wiki/Eccentricity_vector
    '''
    radius = np.linalg.norm(position)
    e = np.linalg.norm(np.cross(velocity, np.cross(position, velocity)) - position / radius)
    return e

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

def import_cpp_data(filename):
    with open(filename, 'r') as file:
        data = np.array(eval(file.read()))
    return data
    

Tmax = 20000
dt = 0.02
nt = int((Tmax - 0) / dt) + 1
NsBH = 10
NsBHMasses = np.ones(NsBH) * 10
SMBHMass = 1e8
r_s = 2 * 4.3e-3 * SMBHMass / 9e10  # 2GM / c^2     units of pc
Nr_s = 1e3      # number of schwarzschild radii to initialise the sim with respect to
lenscale = Nr_s * r_s

times = import_cpp_data("times.txt")
positions = import_cpp_data("positions.txt")
velocities = import_cpp_data("velocities.txt")
nt = len(times)
dt = times[1] - times[0]

real_times = time_convert(times, SMBHMass + sum(NsBHMasses), lenscale)

fig, axes = plt.subplots(figsize=(10, 10), nrows=2, sharex=True, gridspec_kw={'hspace':0})
step = int(nt / 2000) # want 2000 points on our plot
for i in range(1, 11):
    radii = np.array([np.linalg.norm(positions[i, :, j]) for j in range(0, nt, step)])
    vel_mags = np.array([np.linalg.norm(velocities[i, :, j]) for j in range(0, nt, step)])
    semi_majors = - radii / (radii * vel_mags**2 - 2)
    eccentricities = np.array([eccentricity(positions[i, :, j], velocities[i, :, j]) for j in range(0, nt, step)])
    axes[0].plot(real_times[::step], semi_majors * Nr_s)
    axes[1].plot(real_times[::step], eccentricities, lw=0.5)
axes[0].set(yscale='log', ylabel="Semi-Major Axis ($R_s$)")
axes[1].set(yscale='log', xlabel="Time (Myr)", ylabel='Eccentricity')
fig.savefig('NBodyTest.png', dpi=400, bbox_inches='tight')