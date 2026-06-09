%% =========================================================
% TSV RLGC PARAMETER EXTRACTION
% =========================================================

clc;
clear;
close all;

%% =========================================================
% 1. CONSTANTS
% =========================================================

eps0 = 8.854e-12;      % F/m
mu0  = 4*pi*1e-7;      % H/m

%% =========================================================
% 2. INPUT PARAMETERS
% =========================================================

f = 100e9;             % Hz

P_TSV = 10e-6;         % m
d_TSV = 4.94e-6;       % m
r_TSV = 2.47e-6;       % m
h_TSV = 25e-6;         % m

t_ox  = 0.13e-6;       % m
tdep = 0.87e-6         % m

lp    = 12e-6;         % m

theta = pi/2;
J     = 0.10;

T  = 298;              % K
T0 = 273;              % K

%% =========================================================
% 3. MATERIAL PROPERTIES
% =========================================================

sigma_Si  = 10;        % S/m
sigma_TSV = 5.9e7;     % S/m

rho0  = 1.68e-8;       % ohm.m
alpha = 0.039;         % K^-1

mu_r_TSV = 1;

epsr_Si  = 11.9;
epsr_ox  = 4;
epsr_air = 1;

%% =========================================================
% 4. DEFECT-FREE PARAMETERS
% =========================================================

delta_skin = 1/sqrt(pi*f*mu0*mu_r_TSV*sigma_TSV);

Rdc_TSV = rho0*h_TSV/(pi*r_TSV^2);

Rac_TSV = rho0*h_TSV / ...
          (2*pi*r_TSV*delta_skin - pi*delta_skin^2);

RTSV = sqrt(Rdc_TSV^2 + Rac_TSV^2);

LTSV = 0.5 * ...
      ((mu0*mu_r_TSV*h_TSV)/(2*pi)) * ...
      log(P_TSV/r_TSV);

GSi = (pi*sigma_Si*h_TSV) / ...
      acosh(P_TSV/d_TSV);

Cox = 0.5 * ...
      ((pi*eps0*epsr_ox*(h_TSV/2)) / ...
      log((r_TSV+t_ox)/r_TSV));

Cdep = 0.5 * ...
       ((pi*eps0*epsr_Si*(h_TSV/2)) / ...
       log((r_TSV+t_ox+tdep)/(r_TSV+t_ox)));

CSi = (pi*eps0*epsr_Si*h_TSV) / ...
      acosh(P_TSV/d_TSV);

%% =========================================================
% 5. DEFECT PARAMETERS
% =========================================================

Cleak = ...
    (pi*eps0*epsr_ox*(h_TSV/2) ...
    - eps0*epsr_air*theta*lp) ...
    / log((r_TSV+t_ox)/r_TSV);

Gleak = ...
    (2*pi*sigma_Si*h_TSV ...
    - 4*pi^2*J*f*eps0*epsr_air*r_TSV) ...
    / log((r_TSV+t_ox)/r_TSV);

rhoT = rho0*(1 + alpha*(T-T0));

Rsht = rhoT*t_ox / ...
      ((theta/2)*((r_TSV+t_ox)^2 - r_TSV^2));

A = (theta/2)* ...
    ((r_TSV+t_ox)^2 - r_TSV^2);

req = sqrt(A/pi);

Lsht = (mu0*lp/(2*pi))* ...
       (log((2*lp)/req) - 1);

%% =========================================================
% 6. DISPLAY RESULTS
% =========================================================

fprintf('\n');
fprintf('========== DEFECT-FREE PARAMETERS ==========\n');

fprintf('RTSV  = %.6f Ohm\n',RTSV);
fprintf('LTSV  = %.6e H\n',LTSV);
fprintf('GSi   = %.6e S\n',GSi);
fprintf('Cox   = %.6e F\n',Cox);
fprintf('Cdep  = %.6e F\n',Cdep);
fprintf('CSi   = %.6e F\n',CSi);

fprintf('\n');
fprintf('========== DEFECT PARAMETERS ==========\n');

fprintf('Cleak = %.6e F\n',Cleak);
fprintf('Gleak = %.6e S\n',Gleak);
fprintf('Rsht  = %.6f Ohm\n',Rsht);
fprintf('Lsht  = %.6e H\n',Lsht);

fprintf('\n');
fprintf('========== INTERMEDIATE VALUES ==========\n');

fprintf('rhoT  = %.6e Ohm.m\n',rhoT);
fprintf('Area  = %.6e m^2\n',A);
fprintf('req   = %.6e m\n',req);
fprintf('delta = %.6e m\n',delta_skin);