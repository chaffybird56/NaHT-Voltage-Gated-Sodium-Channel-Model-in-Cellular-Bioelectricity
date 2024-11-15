% thresholdcurrent_demo_NaHT.m
% This script determines the threshold current required to trigger an action potential
% with and without the high-threshold Na+ (NaHT) channel.

%% Clear workspace
clc; clearvars; close all;

%% Initialize variables
% max_current and step_size are used in the iterative search for the threshold current.
% Vm0 is the resting membrane potential, and vm is the specific voltage variable for rate calculations.

tstart = 0;               % Initial time in ms
tend = 10;                % Final time in ms
vm=0;                     
Vm0 = -60; 
step_size = 0.1;   
max_current = 50;         % Maximum current to test in μA/cm^2

% The gating particle rates are calculated using Hodgkin-Huxley model equations.
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

% Initial gating particle values
m_0 = alpha_m_0/(alpha_m_0+beta_m_0);
h_0 = alpha_h_0/(alpha_h_0+beta_h_0);
n_0 = alpha_n_0/(alpha_n_0+beta_n_0);
w_0 = alpha_w_0/(alpha_w_0+beta_w_0);
z_0 = alpha_z_0/(alpha_z_0+beta_z_0);

% Initial conditions with NaHT channel
Y0_with_NaHT = [Vm0, m_0, h_0, n_0, w_0, z_0];

% Initial conditions without NaHT channel (without w and z variables)
Y0_without_NaHT = [Vm0, m_0, h_0, n_0, 0, 0];

%% Function for the NaHT channel included
fun_with_NaHT = @(t, Y, I) membranewithNaHT_ODE(t, Y, [I, 1, 0, 0, 0, 40]);  % I, duration 1ms, gNaHT_bar = 40

fun_without_NaHT = @(t, Y, I) membranewithNaHT_ODE(t, Y, [I, 1, 0, 0, 0, 40]);

%% Find threshold current with NaHT channel
current = 0;
threshold_current_with_NaHT = NaN;  % Initialize threshold current as NaN

% Iteratively search for the threshold current required to trigger an
% action potential with the NaHT channel. The loop continues until a
% threshold current is found (threshold_current_with_NaHT is not NaN)
% or the current exceeds the maximum current set for testing.

while isnan(threshold_current_with_NaHT) && current <= max_current
    [t_with_NaHT, Y_with_NaHT] = ode15s(@(t,Y) fun_with_NaHT(t, Y, current), [tstart tend], Y0_with_NaHT);

     % Extract the membrane potential (Vm) from the ODE solution.
    Vm_with_NaHT = Y_with_NaHT(:,1);
    
    % Check if an action potential is triggered by seeing if the membrane potential
    % crosses 0 mV (Vm >= 0) at any point during the simulation.

    if any(Vm_with_NaHT >= 0)  % Assuming an action potential is any Vm crossing 0 mV

    % If an action potential is triggered, record the current as the threshold current.    
        threshold_current_with_NaHT = current;
    else
    
    % If not, increase the current by the predefined step size and try again.
        current = current + step_size;
    end
end

disp(['Threshold current with NaHT channel: ', num2str(threshold_current_with_NaHT), ' μA/cm^2']);

% Determine threshold voltage for with NaHT
peaks_with_NaHT = nonzeros(Vm_with_NaHT.*islocalmax(Vm_with_NaHT)); 
V_threshold_with_NaHT = peaks_with_NaHT(1);  

%% Find threshold current without NaHT channel
current = 0;
threshold_current_without_NaHT = NaN;  % Initialize threshold current as NaN

while isnan(threshold_current_without_NaHT) && current <= max_current
    [t_without_NaHT, Y_without_NaHT] = ode15s(@(t,Y) fun_without_NaHT(t, Y, current), [tstart tend], Y0_without_NaHT);
    Vm_without_NaHT = Y_without_NaHT(:,1);
    
    if any(Vm_without_NaHT >= 0)  % Check if an action potential was triggered
        threshold_current_without_NaHT = current;
    else
        current = current + step_size;
    end
end

disp(['Threshold current without NaHT channel: ', num2str(threshold_current_without_NaHT), ' μA/cm^2']);

%% Find peaks in the membrane potential data to determine the threshold voltage.
% This uses the islocalmax function to identify local maxima, which correspond to action potentials.
peaks_without_NaHT = nonzeros(Vm_without_NaHT.*islocalmax(Vm_without_NaHT));

% The threshold voltage is considered to be the first peak found,
% which represents the voltage at which the first action potential is triggered.
V_threshold_without_NaHT = peaks_without_NaHT(1);  

%% Plot for the case with the NaHT channel included.
figure;
plot(t_with_NaHT, Vm_with_NaHT, 'LineWidth', 2); % Increase line width
title('Threshold Current Response with NaHT Channel');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
hold on;
yline(V_threshold_with_NaHT, '--', 'Threshold with NaHT', 'LineWidth', 2); % Increase line width for threshold line
hold off;

%% Plot for the case without the NaHT channel.
figure;
plot(t_without_NaHT, Vm_without_NaHT, 'LineWidth', 2); % Increase line width
title('Threshold Current Response without NaHT Channel');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
hold on;
yline(V_threshold_without_NaHT, '--', 'Threshold without NaHT', 'LineWidth', 2); % Increase line width for threshold line
hold off;

%% Displaying the threshold potentials for both cases.
disp(['Threshold potential with NaHT: ', num2str(V_threshold_with_NaHT), ' mV']);
disp(['Threshold potential without NaHT: ', num2str(V_threshold_without_NaHT), ' mV']);
