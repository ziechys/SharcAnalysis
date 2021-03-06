# SharcAnalysis

This is a collection of analysis scripts for SHARC, Surface Hopping including Arbitrary Couplings developed by the González group (https://sharc-md.org/).

Some scripts are based on or slight modification of the scripts provided by SHARC (original scripts see: https://github.com/sharc-md/sharc)

## traj_module.py
A python module with functions to read, write and manipulate (Sharc) trajectories.
Based on pymatgen Molecule objects. 

Features: 
- Read xyz trajectory file
- Visualise single Molecules or full trajectories in Jupyter Notebook (using py3Dmol) 
- Operate on molecules and trajectories with symmetry operation
- Molecule matching using Hungarian algorithm
- Combining several trajectory files to combined single trajectory (with weighting) 
- Internal coordinaate analysis of trajectory (or single Molecule)
- FT anaylsis of trajectories (as done in FT_analysis.py) 
- Write trajectory file into separate single Gaussian inputs 
- Read in and plot 1D PES calculations, e.g. from previously calculated points along trajectory 

## FT_analysis.py
Fourier transformation analysis of the relaxation dynamics in SHARC based on structural changes. 
Coded in style with SHARC's original scripts using traj_manip class and geo.py 

Features:
- Define internal coordinates for Fourier Transformation anaylsis
- FT for each trajectory to obtain relevant frequencies in relaxation dynamics based on structural change
- Summed-up FT results over all trajectories
- Creation of coherent wave packet and FT analysis of its structural dynamics
- Various post-processing methods to improve FT
- Plotting of FT results
- pure excited state dynamics and normal mode analysis by detecting (forced) GShop 

Use ```python3 FT_analysis.py``` 

## trajana_nma.py
Normal mode analysis of SHARC trajectories.

Slight modifications of original SHARC script:

- conversion to python3
- fixing of several bugs to make it run
- modification of bar plotting routine
- modification of normal mode analysis (weighted vs. unweighted, coherent vs. total)
- pure excited state dynamics and normal mode analysis by detecting (forced) GShop 
- anharmonic frequency analysis from a Gaussian fchk via Anharmonic Jupyter notebook (integration into main code will come soon)

Use ```python3 trajana_nma.py``` 

## single_traj_analysis.py
Analyse all trajectories of a run individually.
The following analysis is conducted: MCH energies, Diag energies, MCH amplitudes, Diag amplitudes, hopping probabilities, Geometry change (requested through Geo.inp)

Use ```python3 single_traj_analysis.py <path-to-trajectories> <number-of-states> <Geo.inp> <Geo label>``` 

## diagnostics.py
Analyse trajectory runs of SHARC

Slight modifications of original SHARC script:
- Define time window of analysis
- Detect forced GS hop

Use ```python2 diagnostics.py``` 

## setup_traj.py
Setting up trajectories for SHARC run

Slight modifications of original SHARC script:

- fix bug for forced GShop
- include option to specify laserwidth

Use ```python2 setup_traj.py``` 

## corr_laser.py
Script to correct laser hop bug: 
In SHARC resonant laser hops can occur even after the laser intensity reached zero. This script (crudely) solves this problem by altering the SHARC produced (laser.x) laser file and setting the frequency after the laser puls ended to high values and, thus, supressing resonant laser hops.

Use ```python3 corr_laser.py <laser-file> <time after which laser hops should be avoided>``` 

## corr_nan.py
Script to correct 'nan' bug occuring for forced GS hop setup.
Script reads MD output substitutes 'nan' entries in coefficient data with expected coefficient for a GS hop. 

Use ```python3 corr_nan.py <path-to-trajectories>``` 

## end_struc.py

Get all final structures of a MD run in one file (e.g. to view with Molden)

Use ```python3 end_struc.py <path-to-trajectories>``` 

## post_analysis.py (now partially incorporated in diagnostics.py)
Script to analyse SHARC MD simulations with forced GShop.
The script checks if the trajectories behave ordinary befor the GShop and flags them if so or not (DONT_ANALZE flag).

Information output from the script (additional to SHARC's own diagnostics.py): 
- Same gradient in two subsequent time steps during simultion (indicated problem in quantum chemistry calculation)
- When does the GS hop occur
- Is the state changing after the GS hop (it should not)?
- What is T_use from diagnostics.py

Trajectories are groups in:
- no GS hop
- usable trajectory until GS hop, i.e. no criteria set in diagnostics.py or double gradients before GS hop
- not usable trajectories due to T_use > time of GS hop or double gradients before GS hop

Use ```python3 post_analysis.py <path-to-trajectories> <diagnostics.py output>``` 
