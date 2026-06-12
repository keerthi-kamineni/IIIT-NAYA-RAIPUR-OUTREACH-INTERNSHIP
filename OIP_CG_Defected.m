clc;
clear;
close all;

%% Constants
eps0 = 8.854e-12;          % F/m
mu0  = 4*pi*1e-7;          % H/m

%% TSV Parameters
P_TSV  = 10e-6;            % m
d_TSV  = 4.944e-6;         % m
r_TSV  = 2.47e-6;          % m
h_TSV  = 15e-6;            % m
tox    = 0.13e-6;          % m
tdep   = 0.87e-6;          % m

sigma_si   = 10;           % S/m
sigma_TSV  = 5.9e7;        % S/m
sigma_SiO2 = 1e-16;        % S/m

rho_TSV = 1.68e-8;         % ohm-m

er_si  = 11.9;
er_ox  = 4;
er_air = 1;

mur_TSV = 1;

f = 100e9;                 % Hz

theta = pi/2;
Lp = 4e-6;                 % m

J = 0.4;                   % 40%

%% Equation 1 : TSV Resistance

mu_TSV = mu0*mur_TSV;

delta_SD_TSV = 1/sqrt(pi*f*mu_TSV*sigma_TSV);

Rdc_TSV = rho_TSV*h_TSV/(pi*r_TSV^2);

Rac_TSV = rho_TSV*h_TSV/...
          (2*pi*r_TSV*delta_SD_TSV - pi*delta_SD_TSV^2);

R_TSV = sqrt(Rdc_TSV^2 + Rac_TSV^2);

%% Equation 2 : TSV Inductance

L_TSV = 0.5*((mu0*mur_TSV*h_TSV)/(2*pi))*log(P_TSV/r_TSV);

%% Equation 3 : Silicon Conductance

G_si = (pi*sigma_si*(h_TSV-Lp))/acosh(P_TSV/d_TSV);

%% Equation 4 : Oxide Capacitance

Cox = (2*pi*eps0*er_ox*h_TSV)/...
      log((r_TSV+tox)/r_TSV);

%% Equation 5 : Depletion Capacitance

Cdep = (2*pi*eps0*er_si*h_TSV)/...
       log((r_TSV+tox+tdep)/(r_TSV+tox));

%% Equation 6 : Silicon Capacitance

Csi = (pi*eps0*er_si*(h_TSV-Lp))/...
      acosh(P_TSV/d_TSV);

%% Equation 7 : Insulation Capacitance

Cins = (Cdep*Cox)/(Cdep + Cox);

%% Equation 8 : Chao Liu Leakage Capacitance

Cleak = (2*pi*eps0*er_ox*h_TSV ...
        - 2*pi*eps0*er_air*J*r_TSV)/...
        log((r_TSV+tox)/r_TSV);

%% Equation 9 : Chao Liu Leakage Conductance

Gleak = (2*pi*sigma_SiO2*h_TSV ...
        - 4*pi*J*f*eps0*er_air*r_TSV)/...
        log((r_TSV+tox)/r_TSV);

%% Equation 10 : Proposed Short Capacitance

Cair = (theta*eps0*er_air*Lp)/...
       log((r_TSV+tox)/r_TSV);

CsiO2_remain = ((2*pi-theta)*eps0*er_ox*Lp)/...
               log((r_TSV+tox)/r_TSV);

Cshort = CsiO2_remain + Cair;

%% Equation 11 : Modified Oxide Capacitance

Cx = (2*pi*eps0*er_ox*(h_TSV-Lp))/...
     log((r_TSV+tox)/r_TSV);

Cox_new = Cx + Cshort;

%% Equation 12 : Proposed Short Conductance

Gair = (theta*eps0*er_air*Lp*f)/...
       log((r_TSV+tox)/r_TSV);

GsiO2_remain = ((2*pi-theta)*sigma_SiO2*Lp)/...
               log((r_TSV+tox)/r_TSV);

Gshort = GsiO2_remain + Gair;

Gox = Gshort;

%% Display Results

fprintf('R_TSV      = %.6e Ohm\n',R_TSV);
fprintf('L_TSV      = %.6e H\n',L_TSV);
fprintf('G_si       = %.6e S\n',G_si);
fprintf('Cox        = %.6e F\n',Cox);
fprintf('Cdep       = %.6e F\n',Cdep);
fprintf('Csi        = %.6e F\n',Csi);
fprintf('Cins       = %.6e F\n',Cins);
fprintf('Cleak      = %.6e F\n',Cleak);
fprintf('Gleak      = %.6e S\n',Gleak);
fprintf('Cshort     = %.6e F\n',Cshort);
fprintf('Cox_new    = %.6e F\n',Cox_new);
fprintf('Gshort     = %.6e S\n',Gshort);