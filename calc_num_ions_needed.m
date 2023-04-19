% This script calculates the number of ions needed for a desired
% SUPER_ION_MULT value.
clear
close all

%% variables to adjust
super_ion_mult = 300;
HT_CYL_DEBYE = 5;
RAD_CYL_DEBYE = 2;
TEMP_ELC = 58022;
DEN_FAR_PLASMA = 1e14;

%% constants
PERM_FREE_SPACE = 8.84188e-12;
BOLTZMANN = 1.380649e-23;
CHARGE_ELC = 1.602177e-19;

%% derived values
DEBYE = sqrt((PERM_FREE_SPACE * BOLTZMANN * TEMP_ELC)/...
    (DEN_FAR_PLASMA * CHARGE_ELC * CHARGE_ELC));

HT_CYL = HT_CYL_DEBYE * DEBYE;
RAD_CYL = RAD_CYL_DEBYE * DEBYE;

volume = pi * RAD_CYL * RAD_CYL * 2 * HT_CYL;
num_ions = (volume * DEN_FAR_PLASMA) / super_ion_mult;

fprintf(['\n For SUPER_ION_MULT = ' num2str(super_ion_mult)])
fprintf(['\n  min # ions needed = ' num2str(ceil(num_ions)) '\n'])