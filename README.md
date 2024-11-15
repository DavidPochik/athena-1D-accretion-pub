# Accretion problem setup:
 - This version of Athena++ generates parameterized CCSN accretion models.
 - ```jobscript_slurm``` or ```jobscript_local``` run a default accretion model that uses the following:
  - ```Lnu=40``` (in units of $10^{51} \, \mathrm{ergs} \, \mathrm{s}^{-1} \, \mathrm{g}^{-1}$ )
  - ```Mdot=0.7``` (in units of $M_{\odot} \, \mathrm{s}^{-1}$)
  - ```Mach=2.0```
  - ```RPNS=30``` (in units of km)
  - ```enue=12.6156``` (in units of MeV)
  - ```enueb=18.9234``` (in units of MeV) 
 - To see native the Athena++ code and documentation, visit https://github.com/PrincetonUniversity/athena
