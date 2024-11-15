# Accretion problem setup:
	- This version of Athena++ generates parameterized CCSN accretion models.
	- ```jobscript_slurm``` or ```jobscript_local``` run a default accretion model that uses the following:
		1) ```Lnu=40``` (in units of $10^{51}  \mathrm{ergs}  \mathrm{s}^{-1} \mathrm{g}^{-1}$ )
		2) ```Mdot=0.7``` (in units of $M_{\odot} \mathrm{s}^{-1}$)
		3) ```Mach=2.0```
		4) ```RPNS=30``` (in units of km)
		5) ```enue=12.6156``` (in units of MeV)
		6) ```enueb=18.9234``` (in units of MeV) 
		7) ```Tau=2/3```
		8) Boundary conditions described in Pochik and Thompson 2025 (in progress)
	- Variables in the ```athinput``` files already include the default boundary values and conditions listed above.
	- The ```accretion_12panel.py``` plot script produces a plot of the accretion model at the final timestep.
	- To see native the Athena++ code and documentation, visit https://github.com/PrincetonUniversity/athena
