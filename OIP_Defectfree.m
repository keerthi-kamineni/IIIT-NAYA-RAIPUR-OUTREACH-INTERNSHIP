 clc;
clear;
close all;

%% =========================================================
% CONSTANTS
% ==========================================================

mu0  = 4*pi*1e-7;        % H/m
eps0 = 8.854e-12;        % F/m

%% =========================================================
% INPUT PARAMETERS
% ==========================================================

P_TSV = 10e-6;           % m
d_TSV = 4.94e-6;         % m
r_TSV = 2.47e-6;         % m
h_TSV = 15e-6;           % m

tox  = 0.13e-6;          % m
tdep = 0.83e-6;          % m
lsub = 3.06e-6;          % m

f = 100e9;               % Hz

%% =========================================================
% MATERIAL PROPERTIES
% ==========================================================

sigma_Si  = 10;          % S/m
sigma_TSV = 5.9e7;       % S/m

rho_TSV = 1.68e-8;       % Ohm-m

epsr_Si = 11.9;
epsr_ox = 4.0;

mur_TSV = 1;

%% =========================================================
% EQUATION 1 : SKIN DEPTH
% ==========================================================

mu_TSV = mu0 * mur_TSV;

delta_SD = 1 / sqrt(pi * f * mu_TSV * sigma_TSV);

%% =========================================================
% EQUATION 2 : TSV RESISTANCE
% ==========================================================

Rdc_TSV = (rho_TSV * h_TSV) / (pi * r_TSV^2);

Rac_TSV = (rho_TSV * h_TSV) / ...
          (2*pi*r_TSV*delta_SD - pi*delta_SD^2);

R_TSV = sqrt(Rdc_TSV^2 + Rac_TSV^2);

%% =========================================================
% EQUATION 3 : TSV INDUCTANCE
% ==========================================================

L_TSV = 0.5 * ((mu0*mur_TSV*h_TSV)/(2*pi)) ...
        * log(P_TSV/r_TSV);

%% =========================================================
% EQUATION 4 : SILICON CONDUCTANCE
% ==========================================================

G_Si = (pi*sigma_Si*h_TSV) / ...
       acosh(P_TSV/d_TSV);

%% =========================================================
% EQUATION 5 : OXIDE CAPACITANCE
% ==========================================================

Cox = (2*pi*eps0*epsr_ox*h_TSV) / ...
      log((r_TSV + tox)/r_TSV);

%% =========================================================
% EQUATION 6 : DEPLETION CAPACITANCE
% ==========================================================

Cdep = (2*pi*eps0*epsr_Si*h_TSV) / ...
       log((r_TSV + tox + tdep)/(r_TSV + tox));

%% =========================================================
% EQUATION 7 : SILICON CAPACITANCE
% ==========================================================

CSi = (pi*eps0*epsr_Si*h_TSV) / ...
      acosh(P_TSV/d_TSV);

%% =========================================================
% EQUATION 8 : INSULATION CAPACITANCE
% ==========================================================

Cins = (Cdep * Cox)/(Cdep + Cox);

%% =========================================================
% DISPLAY RESULTS
% ==========================================================

fprintf('\n========================================\n');
fprintf('      DEFECT-FREE TSV PARAMETERS\n');
fprintf('========================================\n\n');

fprintf('Skin Depth (delta_SD) = %.6e m\n\n', delta_SD);

fprintf('Rdc_TSV = %.6e Ohm\n', Rdc_TSV);
fprintf('Rac_TSV = %.6e Ohm\n', Rac_TSV);
fprintf('R_TSV   = %.6e Ohm\n\n', R_TSV);

fprintf('L_TSV   = %.6e H\n\n', L_TSV);

fprintf('G_Si    = %.6e S\n\n', G_Si);

fprintf('Cox     = %.6e F\n', Cox);
fprintf('Cdep    = %.6e F\n', Cdep);
fprintf('CSi     = %.6e F\n', CSi);
fprintf('Cins    = %.6e F\n', Cins);

fprintf('\n========================================\n');