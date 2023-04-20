# DRIAD_analysis

This repository contains MATLAB code for visualizing the output data from DRIAD. The code is organized into separate files for each type of plot, with a main script that calls the individual plot files to generate all the plots in one run. In total, 14 plots should be made, but each script can be edited to only create plots for what is needed. The main code has the ability to read in multiple input datasets and will save each plot in its dataset folder.

## Requirements

The code is written in MATLAB and requires a MATLAB installation to run. Additionally, the code requires the following MATLAB toolboxes:

- Signal Processing Toolbox
- Optimization Toolbox
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox
- Usage

To use the code, simply clone the repository and open MATLAB. Run the "plots.m" script to generate all the plots for a given simulation dataset. 

## File structure

- plots.m: main script that calls individual plot files to generate all plots.
- import_debug.m: read in data from the debug file.
- plot_potentialscontour.m: plot potential map.
- plot_densitycontour.m: plot density map.
- plot_den_vel_acc.m: plot ion density, velocity, and acceleration.
- dust_charge_avg.m: calculate and plot final dust charge for dust-specific cases.
- charge_time.m: plot dust charge as time passes for dust-specific cases.
- forces.m: plot Fid and Fdd for each dust grain for dust-specific cases.
- import_without_A.m: read in data from the debug file.
- calc_num_ions_needed.m: calculates how many ions needed for a desired SUPER_ION_MULT, given HT_CYL and RAD_CYL.
- plot_3Dionpaths.m: plots the entrance point, trajectory, and exit points of ions in the simulation cylinder.
- drift_velocity.m: calculates the input Mach for a given electric field.

## File input
_ion_on_dust_acc.txt: contains "accDustIon" in x-, y-, and z-directions, followed by "momIonDust" in x-, y-, and z-directions (for each dust grain at each time step).

_debug.txt: contains GPU device properties, constants, user-defined parameters, derived parameters, as well as other values.

_dust-pos.txt: contains the final dust positions.

_params.txt: contains all input parameters.

_ion-den.txt: contains the x and z grid positions in two columns, followed by the ion density in one column and ion potential at each grid position in the other.

_outside_potential.txt: contains the x and z grid positions in two columns, followed by the outside potential at each grid position.
