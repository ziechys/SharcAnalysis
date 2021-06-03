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

## setup_traj.py
Setting up trajectories for SHARC run

Slight modifications of original SHARC script:

- fix big for forced GShop
- include option to specify laserwidth
