% This script calculates the drift velocity for ions entering the
% simulation, given a constant, linear, or quadratic electric field.
% From Nanbu and Kitatani, J. Phys D, 28,324, 1995
% Notes: Lab Notebook #3 page 121

clear
%% SI units
k = 1.38e-23; % Boltzmann constant in J/K
kb_eV = 8.6173303e-5;% Boltzmann constant eV/K

% For argon:
mi = 6.68e-26; % mass in kg for argon
Ti = 290; % K
Te = 58022; % K
alpha_ion = 1.642e-30; %Argon polarizability in m^3

% For neon:
%mi = 3.38e-26; % mass in kg for neon
%Ti = 465.9;% K
%Te = 29011;
%Te = 40615;
%Te = 50107;
%Te = 5/kb_eV; % temperature of electrons eV -> K
%alpha_ion = 0.392e-30; %Neon polarizability in m^3

A = 2.6*(1/1.602e-19)^(1/4);
q = 1.602e-19; %charge in statcoulombs
epsilon0 = 8.85e-12; %SI units
fourpiepsilon0 = 4*pi*epsilon0; %

alpha_ion = 1.642e-30; %Argon polarizability in m^3
a_ion = alpha_ion*q^2/(2*fourpiepsilon0);
bex_ion = A * (4*a_ion)^(1/4);
%bex_neon = 154e-12;

K = 8/3/q * sqrt(mi/pi)*pi*(bex_ion)^2;
z = 0.00162;
E_o = -5817;
alpha = 2.666e+61; 
beta = 0;
E = E_o + alpha*z + beta*z*z;
kb = 1.38e-23; %SI units
P = 19.99; %pressure in Pa
N = P/kb/Ti;% m^-3 
  
vd=(((3*k*Ti/mi)^2 + 3/K^2/mi * (E/N)^2)^(1/2) - 3*k*Ti/mi)^(1/2)

Mach_thermal = vd/sqrt(k*Ti/mi);
Mach_sound_speed = vd/sqrt(k*Te/mi) %assuming 5eV electrons