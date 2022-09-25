# DaMaSCUS-EarthCapture

This is a code to compute the fraction of dark matter that is captured inside the Earth through spin-independent or spin-dependent scattering.

## Instructions for DaMaSCUS-EarthCapture

The code is developed based on [DaMaSCUS](https://github.com/temken/DaMaSCUS). The main code is located in the src/simulation and src/simulation_SD folders.

Requirements:
- openmpi
- eigen3
- boost

To run the code for spin-independent interactions,
```
mpirun -n NCPU ./DaMaSCUS-Simulator config.cfg mDM xs Nsim 
```
where the different parameters are:
- NCPU: number of processors in computation.
- config.cfg: the configuration file.
- mass: the mass of dark matter in MeV.
- xs: dark matter-nucleon scattering cross section in pb.
- Nsim: total number of dark matter particles to be simulated.

If mass, xs, or Nsim is not specified, their values will be set by the corresponding entries in config.cfg.
If mpi error is encounter, one can try
```
 export TMPDIR=/tmp before the run
 ```
 before the run.
 
 To run the code for spin-dependent interactions,
 ```
mpirun -n NCPU ./DaMaSCUS-SD config.cfg mDM xs Nsim 
```
It is also possible to set proton-only or neutron-only scattering in src/simulation_SD/PREM.cpp with
```
#define ap 0
```
or
```
#define an 0
```
To get rid of the solar gravitational acceleration of dark matter velocity, simply set 
```
double vsunEarth=0;
```
in IC_Generator.cpp.

## Using the code:

Feel free to use, modify or distribute the code. If you use the code in your publication, please cite the paper [...](...)

In addition, please cite the original DaMaSCUS code [https://github.com/temken/DaMaSCUS](https://github.com/temken/DaMaSCUS)

For questions or suggestions, please contact Ningqiang Song [ningqiang.song@liverpool.ac.uk](ningqiang.song@liverpool.ac.uk)