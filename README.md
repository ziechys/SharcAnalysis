# SharcAnalysis

This is a collection of analysis scripts for SHARC, Surface Hopping including Arbitrary Couplings developed by the González group (https://sharc-md.org/).

Some scripts are based on or slight modification of the scripts provided by SHARC

## trajana_nma.py
Normal mode analysis of SHARC trajectories.

Slight modifications of original SHARC script:

- conversion to python3
- fixing of several bugs to make it run
- modification of bar plotting routine
- modification of normal mode analysis (weighted vs. unweighted, coherent vs. total)
- pure excited state dynamics and normal mode analysis by detecting GShop (forced)

Use ```python3 trajana_nma.py``` 

## single_traj_analysis.py
Analyse all trajectories of a run individually.
The following analysis is conducted: MCH energies, Diag energies, MCH amplitudes, Diag amplitudes, hopping probabilities, Geometry change (requested through Geo.inp)

Use ```python3 single_traj_analysis.py <path-to-trajectories> <number-of-states> <Geo.inp>``` 

## setup_traj.py
Setting up trajectories for SHARC run

Slight modifications of original SHARC script:

- fix bug for forced GShop
- include option to specify laserwidth

Use ```python2 setup_traj.py``` 

## corr_nan.py
Script to correct 'nan' bug occuring for forced GS hop setup.
Script reads MD output substitutes 'nan' entries in coefficient data with expected coefficient for a GS hop. 

Use ```python3 corr_nan.py <path-to-trajectories>``` 

## end_struc.py

Get all final structures of a MD run in one file (e.g. to view with Molden)

Use ```python3 end_struc.py <path-to-trajectories>``` 
