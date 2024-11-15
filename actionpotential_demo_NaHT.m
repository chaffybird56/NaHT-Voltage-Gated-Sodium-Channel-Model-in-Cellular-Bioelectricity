%% Clear workspace
% This ensures a clean starting environment by clearing variables, closing figures, etc.
clc; clearvars; close all;

%% Initialize variables
% Setting up the simulation parameters and initial conditions for the neuronal model.
tstart = 0;     % Start time of the simulation in milliseconds (ms)
tend = 100;     % End time of the simulation in ms
Vm0 = -60;      % Resting membrane potential in millivolts (mV)
vm = 0;         % Membrane potential relative to threshold (Vh) for rate calculations

% Current injection parameters for simulations with and without the high-threshold Na+ channels (NaHT)
I_parameters_with_NaHT = [20, tend, 0, 0, 0, 40]; % Current step with g_NaHT set to 40 mS/cm^2
I_parameters_without_NaHT = [20, tend, 0, 0, 0, 0]; % Current step with g_NaHT set to 0, channel blocked

%% Calculate gating particle transition rates based on initial conditions
% These parameters define the rate at which ion channels change their states
% from closed to open and vice versa, under the initial membrane potential conditions.
alpha_m_0 = 0.374*(vm-25.41)/(1-exp((25.41-vm)/6.06));
beta_m_0 = 0.795*(21-vm)/(1-exp((vm-21)/9.41));
alpha_h_0 = -0.110*(27.74+vm)/(1-exp((27.74+vm)/9.06));
beta_h_0 = 4.514/(1+exp((56-vm)/12.5));
alpha_n_0 = 0.0516*(vm-35)/(1-exp((35-vm)/10));
beta_n_0 = 0.129*(35-vm)/(1-exp((vm-35)/10));
alpha_w_0 = 0.0936 * (vm - 55.41) / (1 - exp((55.41 - vm) / 6.06));
beta_w_0 = 0.199 * (51 - vm) / (1 - exp((vm - 51) / 9.41));
alpha_z_0 = -0.055 * (27.74 + vm) / (1 - exp((vm + 27.74) / 9.06));
beta_z_0 = 2.257 / (1 + exp((56 - vm) / 12.5));

%% Compute initial gating particle values
% These values represent the initial state of each gating particle.
% They are calculated from the alpha and beta transition rates.
m_0 = alpha_m_0/(alpha_m_0+beta_m_0);
h_0 = alpha_h_0/(alpha_h_0+beta_h_0);
n_0 = alpha_n_0/(alpha_n_0+beta_n_0);
w_0 = alpha_w_0/(alpha_w_0+beta_w_0);
z_0 = alpha_z_0/(alpha_z_0+beta_z_0);

% Combine all initial states into a vector for use in ODE solving.
Y0 = [Vm0, m_0, h_0, n_0, w_0, z_0];

%% Set ODE solver options
% These options help control the precision and step size of the ODE solver.
options = odeset('MaxStep', 1e-3, 'RelTol', 1e-3);

%% Solve ODE for the case WITH NaHT channel
[t_with_NaHT, Y_with_NaHT] = ode15s(@(t,Y) membranewithNaHT_ODE(t, Y, I_parameters_with_NaHT), [tstart tend], Y0, options);
Vm_with_NaHT = Y_with_NaHT(:,1); % Membrane potential with NaHT channel

%% Solve ODE for the case WITHOUT NaHT channel
[t_noNaHT, Y_noNaHT] = ode15s(@(t,Y) membranewithNaHT_ODE(t, Y, I_parameters_without_NaHT), [tstart tend], Y0, options);
Vm_noNaHT = Y_noNaHT(:,1); % Membrane potential without NaHT channel

%% Analysis of spike rates
% Using 'islocalmax' to detect local maxima in the membrane potential,
% which correspond to action potentials (spikes).

% Analysis WITH NaHT channel
peaks_with_NaHT = nonzeros(t_with_NaHT.*islocalmax(Vm_with_NaHT)); 
N_with_NaHT = length(peaks_with_NaHT);
av_period_with_NaHT = mean(diff(peaks_with_NaHT)); 
spike_rate_with_NaHT = 1 / av_period_with_NaHT * 1000; % Convert to spikes per second (Hz)
disp("Spike rate WITH NaHT: " + string(spike_rate_with_NaHT) + " spikes/s");

% Analysis WITHOUT NaHT channel
peaks_noNaHT = nonzeros(t_noNaHT.*islocalmax(Vm_noNaHT)); 
N_noNaHT = length(peaks_noNaHT);
av_period_noNaHT = mean(diff(peaks_noNaHT)); 
spike_rate_noNaHT = 1 / av_period_noNaHT * 1000; % Convert to spikes per second (Hz)
disp("Spike rate WITHOUT NaHT: " + string(spike_rate_noNaHT) + " spikes/s");

%% Plotting results
% Plot for WITH NaHT channel
figure;
plot(t_with_NaHT, Vm_with_NaHT, 'LineWidth', 2); % Setting linewidth to 2
title('Membrane Potential with NaHT Channel');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');

% Plot for WITHOUT NaHT channel
figure;
plot(t_noNaHT, Vm_noNaHT, 'LineWidth', 2); % Setting linewidth to 2
title('Membrane Potential without NaHT Channel');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
