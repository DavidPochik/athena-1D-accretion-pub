# Accretion problem setup:
 - This version of Athena++ generates parameterized CCSN accretion models.
 - ```jobscript_parallel_main``` or ```jobscript_serial_main``` run a default accretion model that uses the following:
  1) ```Lnu=40``` (in units of $10^{51}  \mathrm{ergs}  \mathrm{s}^{-1} \mathrm{g}^{-1}$ )
  2) ```Mdot=0.7``` (in units of $M_{\odot} \mathrm{s}^{-1}$)
  3) ```Mach=2.0```
  4) ```RPNS=30``` (in units of km)
  5) ```enue=12.6156``` (in units of MeV)
  6) ```enueb=18.9234``` (in units of MeV) 
  7) ```Tau=2/3```
  8) Boundary conditions described in Pochik and Thompson 2025 (https://arxiv.org/abs/2411.16857)
 - Variables in the ```athinput``` files already include the default boundary values and conditions listed above.
 - The ```accretion_12panel.py``` plot script produces a plot of the accretion model at the final timestep.
 - ```jobscript_parallel_thermal``` runs the code used for the accretion model analysis in Section 3.6 of Pochik and Thompson 2025 using the following:
  1) ```Lnu=10``` (in units of $10^{51}  \mathrm{ergs}  \mathrm{s}^{-1} \mathrm{g}^{-1}$ )
  2) ```Mdot=0.3``` (in units of $M_{\odot} \mathrm{s}^{-1}$)
  3) ```S_out=3.``` (in units of $k_{\mathrm{B}}$ baryon$^{-1}$)
  4) ```vcoef=0.8``` (fraction of free fall)
  5) ```Ye_out=26/56```
  6) ```Abar_out=56.```
  7) ```RPNS=30``` (in units of km)
  8) ```MPNS=1.4``` (in units of $M_{\odot}$)
  9) ```enue=12.6156``` (in units of MeV)
  10) ```enueb=18.9234``` (in units of MeV) 
  11) ```Tau=2/3``` 
- To see native the Athena++ code and documentation, visit https://github.com/PrincetonUniversity/athena

