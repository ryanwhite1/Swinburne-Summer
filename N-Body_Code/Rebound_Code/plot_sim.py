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

# Tmax = 20000
# dt = 0.02
# nt = int((Tmax - 0) / dt) + 1
# NsBH = 10
# NsBHMasses = np.ones(NsBH) * 10
# SMBHMass = 1e8
# r_s = 2 * 4.3e-3 * SMBHMass / 9e10  # 2GM / c^2     units of pc
# Nr_s = 1e3      # number of schwarzschild radii to initialise the sim with respect to
# lenscale = Nr_s * r_s

# rawdata = np.genfromtxt('orbits.txt')
# r, c = rawdata.shape
# rawdata = rawdata.reshape((10, r//10, c), order='F')

# times = rawdata[0, :, 0]
# # # positions = import_cpp_data("positions.txt")
# # # velocities = import_cpp_data("velocities.txt")
# nt = len(times)
# dt = times[1] - times[0]

# real_times = time_convert(times, SMBHMass + sum(NsBHMasses), lenscale)

# fig, axes = plt.subplots(figsize=(10, 10), nrows=2, sharex=True, gridspec_kw={'hspace':0})
# for i in range(10):
#     semi_majors = rawdata[i, :, 1]
#     eccentricities = rawdata[i, :, 2]
#     axes[0].plot(real_times, semi_majors * Nr_s)
#     axes[1].plot(real_times, eccentricities, lw=0.5)
# axes[0].set(yscale='log', ylabel="Semi-Major Axis ($R_s$)")
# axes[1].set(yscale='log', xlabel="Time (Myr)", ylabel='Eccentricity')
# fig.savefig('NBodyTest.png', dpi=400, bbox_inches='tight')
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.cm import ScalarMappable

SMBHMass = 1e8
r_s = 2 * 4.3e-3 * SMBHMass / 9e10  # 2GM / c^2     units of pc
Nr_s = 1e3      # number of schwarzschild radii to initialise the sim with respect to
lenscale = Nr_s * r_s
rawdata = np.genfromtxt('orbits.txt')
N = int(rawdata[:, 1].max())
fig, axes = plt.subplots(figsize=(10, 10), nrows=2, sharex=True, gridspec_kw={'hspace':0, 'height_ratios':[1.5, 1]})

# rawdata[:, 5] *= SMBHMass
# norm = plt.Normalize(rawdata[:, 5].min(), rawdata[:, 5].max())
# cmap = 'viridis'
fig2, ax2 = plt.subplots()

for i in range(1, N+1):
    particle_data = rawdata[np.where(rawdata[:, 1] == i)[0], :]
    semi_majors = particle_data[:, 2] * Nr_s
    eccentricities = particle_data[:, 3]
    inclinations = particle_data[:, 4]  # convert from rad to degrees
    real_times = time_convert(particle_data[:, 0], SMBHMass, lenscale)
    line, = axes[0].plot(real_times, semi_majors, lw=0.5)
    axes[1].plot(real_times, eccentricities, lw=0.5)
    ax2.plot(real_times, inclinations, lw=0.5)
    for j in [0, 1]:
        axes[j].axvline(real_times[-1], ls='--', lw=0.5, c='k')
    ax2.axvline(real_times[-1], ls='--', lw=0.5, c='k')
        
    # points = np.array([real_times, semi_majors * Nr_s]).T.reshape(-1, 1, 2)
    # segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # lc = LineCollection(segments, cmap=cmap, norm=norm)
    # lc.set_array(particle_data[:, 5])
    # lc.set_linewidth(1)
    # line = axes[0].add_collection(lc)
    
    ### now to overlay the mass changes
    masses = particle_data[:, 5] * SMBHMass
    axes[0].text(real_times[0], semi_majors[0], f"{masses[0]:.1f}", c=line.get_color(), fontsize=4, ha='right')
    unique_masses = np.unique(masses)
    if len(unique_masses) > 1:
        for j in range(1, len(unique_masses)):
            index = np.argwhere(masses == unique_masses[j]).flatten()[0]
            mass_text = f"{masses[index - 1]:.1f} -> {masses[index]:.1f}"
            axes[0].text(real_times[index], semi_majors[index], mass_text, c=line.get_color(), fontsize=4, ha='center')
    
    
axes[0].set(yscale='log', ylabel="Semi-Major Axis ($R_s$)")
axes[1].set(yscale='log', xlabel="Time (Myr)", ylabel='Eccentricity')
# fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes[0], label='Mass ($M_\odot$)', 
#              location='top', orientation='horizontal', aspect=50, pad=0)
fig.savefig('NBody_a-e_Plot.png', dpi=500, bbox_inches='tight')
ax2.set(xlabel="Time (Myr)", ylabel='Inclination (degrees)', yscale='log')
fig2.savefig('NBody_inclination_Plot.png', dpi=500, bbox_inches='tight')


### now to plot the binary mass and binary mass ratio plot
merger_data = np.genfromtxt('Mergers.txt')
mergers = merger_data.shape[0]

binary_masses = np.array([merger_data[i, 1] + merger_data[i, 2] for i in range(mergers)])
binary_m_ratio = np.array([min(merger_data[i, 1], merger_data[i, 2]) / max(merger_data[i, 1], merger_data[i, 2]) for i in range(mergers)])
effective_spins = merger_data[:, 4]
remnant_spins = merger_data[:, 5]

fig, ax = plt.subplots()
n_bins = 10
xbins = np.logspace(np.log10(min(binary_masses)), np.log10(max(binary_masses)), n_bins)
ybins = np.logspace(np.log10(min(binary_m_ratio)), np.log10(max(binary_m_ratio)), n_bins)
_, _, _, cbar = ax.hist2d(binary_masses, binary_m_ratio, cmin=1, bins=[xbins, ybins])
ax.set(xscale='log', yscale='log', xlabel='Binary Mass, $m_1 + m_2$ ($M_\odot$)', ylabel='Mass Ratio $q$ ($m_1 / m_2$)')
fig.colorbar(cbar, label='Counts')
fig.savefig('Q_vs_BinaryMass.png', dpi=400, bbox_inches='tight')


fig, ax = plt.subplots()
ax.hist(remnant_spins, bins=20, ec='k')
ax.set(xlabel="Dimensionless Spin Parameter", ylabel='Frequency')
fig.savefig('BH_Spins.png', dpi=400, bbox_inches='tight')


def callister(x, y):
    # q-vs-chi_eff probability density function from callister (2021) - Who Ordered That? Unequal-mass Binary Black Hole Mergers Have Larger Effective Spins
    return y**1.08 * np.exp(-(x - (0.19 - 0.46*(y - 0.5)))**2 / (2 * (10**(-1.06 - 0.83*(y - 0.5)))**2))
fig, ax = plt.subplots()
x = np.linspace(-1, 1, 100)
y = np.linspace(0, 1, 100)
z = np.ndarray((len(x), len(y)))
for i in range(len(x)):
    for j in range(len(y)):
        z[i, j] = callister(x[i], y[j])
ax.contourf(x, y, z.T, cmap='binary', alpha=0.3, levels=[0.01, 0.05, 0.16, 0.5, 0.84, 1], antialiased=True)
cbar = ax.scatter(effective_spins, binary_m_ratio, c=merger_data[:, 8])
ax.set(xlabel='$\chi_{eff}$', ylabel='Mass Ratio $q$ ($m_1 / m_2$)', xlim=[-1, 1], ylim=[0, 1])
# ax.legend()
fig.colorbar(cbar, label='BH Generation')
fig.savefig('BH_Effective_Spins.png', dpi=400, bbox_inches='tight')

plt.close('all')
