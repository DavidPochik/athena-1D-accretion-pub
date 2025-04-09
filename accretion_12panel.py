import matplotlib
from matplotlib.ticker import FixedLocator
import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.insert(0, str(sys.argv[2]))
import athena_read
try:
	import helmeos
except ImportError:
	try:
		from .. import helmeos
	except ImportError:
		import sys
		sys.path.append(os.path.dirname(os.path.abspath(__file__)))
		import helmeos

DIR_name1 = str(sys.argv[1])
rmin      = 3.0e6 # minimum radius for plots, in cm
rmax      = 1.5e7 # maximum radius for plots, in cm
Mdot      = 0.7
Mach      = 2.0
Lnu       = 40
RPNS      = 30

def dataread(DIR_name):
	n1     = DIR_name+"/accretion.prim.00100.athdf"
	n2     = DIR_name+"/accretion.uov.00100.athdf"
	data1  = athena_read.athdf(n1)
	data2  = athena_read.athdf(n2)

	r     = data1['x1v']
	v     = (data1['vel1'])[0][0]
	p     = (data1['press'])[0][0]
	rho   = (data1['rho'])[0][0]
	T     = (data2['dt1'])[0][0]
	qdot  = (data2['dt2'])[0][0]
	Ye    = (data1['r0'])[0][0]
	cs    = (data2['dt3'])[0][0]
	nR    = np.size(r)

	abar  = 1.0
	G     = 6.6743*pow(10,-8)          # cm^3 g^-1 s^-2
	kb    = 8.61733326e-11 # MeV/K
	M     = 1.4*2*pow(10,33)           # g
	mu    = G*M
	nR    = np.size(r)
	Bern1 = np.zeros(nR)
	Vesc  = np.zeros(nR)
	ante  = np.zeros(nR)

	fn          = os.path.join(os.path.dirname(os.path.abspath(__file__)), str(sys.argv[3]))
	by_filename = helmeos.HelmTable(fn=fn, temp_n=201, dens_n=541)
	EOSData1    = by_filename.eos_DT(rho, T, abar, Ye)
	etot1       = EOSData1['etot'] # erg/g

	for j in range(0,nR):
		Bern1[j] = 0.5 * v[j]**2  + etot1[j] + p[j]  / rho[j]  - mu / r[j]
		Vesc[j]  = np.sqrt(2.0 * G * M / (r[j]))
		ante[j]  = cs[j]**2 / Vesc[j]**2

	return [r, rho, v, T*kb, Ye, qdot, ante, etot1, Bern1, p]

[r1, rho1, v1, T1, Ye1, qdot1, ante1, etot1, Bern1, P1] = dataread(DIR_name1) # Lcore=2

plt.clf()
axs, fig = plt.subplots(3,3, figsize=(9.5,9.5))
namemulti = 'accretion_multiplot'
yfs       = 16
xfs       = 16
fs_tick   = 14
nticksy   = 7
nticksy2  = 5

# Density
fig[0,0].plot(r1/1.0e5,rho1,linewidth=1.25,linestyle='solid',color='black')
fig[0,0].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[0,0].set_yscale('log')
fig[0,0].set_ylabel(r'$\rho$ [g cm$^{-3}$]', fontsize=yfs)
maj_loc00y = ticker.LogLocator(numticks=nticksy)
min_loc00y = ticker.LogLocator(subs='all',numticks=nticksy)
fig[0,0].yaxis.set_major_locator(maj_loc00y)
fig[0,0].yaxis.set_minor_locator(min_loc00y)
fig[0,0].tick_params(axis='x', length=1,labelsize=fs_tick)
fig[0,0].tick_params(axis='y', length=1, labelsize=fs_tick)

# Velocity
fig[0,1].plot(r1/1.0e5,-1.0 * v1,linewidth=1.25,linestyle='solid',color='black')
fig[0,1].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[0,1].set_ylabel(r'$-v_{r}$ [cm s$^{-1}$]', fontsize=yfs)
fig[0,1].tick_params(axis='x', labelsize=fs_tick)
fig[0,1].tick_params(axis='y', labelsize=fs_tick)

# Temperature
fig[0,2].plot(r1/1.0e5,T1,linewidth=1.25,linestyle='solid',color='black')
fig[0,2].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[0,2].set_yscale('log')
fig[0,2].set_xscale('linear')
fig[0,2].set_ylabel(r'$T$ [MeV]', fontsize=yfs)
fig[0,2].tick_params(axis='x', labelsize=fs_tick)
fig[0,2].tick_params(axis='y', labelsize=fs_tick)

# Electron Fraction
fig[1,0].plot(r1/1.0e5,Ye1,linewidth=1.25,linestyle='solid',color='black')
fig[1,0].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[1,0].set_ylabel(r'$Y_{e}$', fontsize=yfs)
fig[1,0].tick_params(axis='x', labelsize=fs_tick)
fig[1,0].tick_params(axis='y', labelsize=fs_tick)

# Pressure
fig[1,1].plot(r1/1.0e5,P1, linewidth=1.25,linestyle='solid',color='black')
fig[1,1].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[1,1].set_yscale('log')
fig[1,1].set_ylabel(r'$P$ [$\mathrm{ergs} \, \mathrm{cm}^{-3}$]', fontsize=yfs)
maj_loc11y = ticker.LogLocator(numticks=nticksy2)
min_loc11y = ticker.LogLocator(subs='all',numticks=nticksy2)
fig[1,1].yaxis.set_major_locator(maj_loc11y)
fig[1,1].yaxis.set_minor_locator(min_loc11y)
fig[1,1].tick_params(axis='x', labelsize=fs_tick)
fig[1,1].tick_params(axis='y', labelsize=fs_tick)

# Qdot
fig[1,2].plot(r1/1.0e5,qdot1/1.0e21,linewidth=1.25,linestyle='solid',color='black')
fig[1,2].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[1,2].set_ylabel(r'$\dot{q}$ [$10^{21}$ erg s$^{-1}$ g$^{-1}$]', fontsize=yfs)
fig[1,2].tick_params(axis='x', labelsize=fs_tick)
fig[1,2].tick_params(axis='y', labelsize=fs_tick)

# Antesonic ratio
fig[2,0].plot(r1/1.0e5,ante1,linewidth=1.25,linestyle='solid',color='black')
fig[2,0].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[2,0].set_ylabel(r'$c_{s}^{2}/v_{\mathrm{esc}}^{2}$', fontsize=yfs)
fig[2,0].set_xlabel(r'$R$ [km]', fontsize=xfs)
fig[2,0].tick_params(axis='x', labelsize=fs_tick)
fig[2,0].tick_params(axis='y', labelsize=fs_tick)

# Internal energy
fig[2,1].plot(r1/1.0e5,etot1,linewidth=1.25,linestyle='solid',color='black')
fig[2,1].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[2,1].set_ylabel(r'$\epsilon$ [erg g$^{-1}$]', fontsize=yfs)
fig[2,1].set_xlabel(r'$R$ [km]', fontsize=xfs)
fig[2,1].tick_params(axis='x', labelsize=fs_tick)
fig[2,1].tick_params(axis='y', labelsize=fs_tick)

# Bernoulli integral
fig[2,2].plot(r1/1.0e5,Bern1/1.0e19,linewidth=1.25,linestyle='solid',color='black')
fig[2,2].set_xlim([rmin/1.0e5, rmax/1.0e5])
fig[2,2].set_ylabel(r'$\mathcal{B}$ [$10^{19}$ erg g$^{-1}$]', fontsize=yfs)
fig[2,2].set_xlabel(r'$R$ [km]', fontsize=xfs)
fig[2,2].tick_params(axis='x', labelsize=fs_tick)
fig[2,2].tick_params(axis='y', labelsize=fs_tick)

plt.tight_layout(pad=4.0)
plt.savefig(namemulti+'.png', bbox_inches='tight')
