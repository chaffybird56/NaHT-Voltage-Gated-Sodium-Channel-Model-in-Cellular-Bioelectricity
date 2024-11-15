% refractory_demo_NaHT.m
% This script demonstrates the relative refractory period properties of the NaHT model
% and the blocked-channel model for a pair of 1-ms long current pulses.

%% Clear workspace
clc; clearvars; close all;

%% Initialize variables
tstart = 0;    % Start time in ms
tend = 50;     % End time in ms
Vm0 = -60;     % Resting membrane potential in mV
gap_blocked = [24.5, 25];  % Inter-pulse gap in ms for blocked NaHT channel
gap_unblocked = [26.5, 27]; % Inter-pulse gap in ms for unblocked NaHT channel
vm = 0; 

% Parameters to find refractory period
initial_gap = 5;  % An initial guess after the absolute refractory period
max_gap = 100;     % Maximum gap to test
gap_increment = 0.1;  % Increment gap in 0.1 ms steps

% Initialize variables for tracking refractory periods
refractory_period_with_NaHT = NaN;
refractory_period_without_NaHT = NaN;
AP_triggered_with_NaHT = false;
AP_triggered_without_NaHT = false;

%% Threshold currents & potentials obtained from previous analysis (q3)
I_threshold_blocked = 23.7; 
I_threshold_unblocked = 22.4; 

V_threshold_with_NaHT = -43.5934;
V_threshold_without_NaHT = -43.2337;

%% Initial conditions 

% Initial gating particle transition rates
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

% Initial conditions without NaHT channel (no w and z variables)
Y0_without_NaHT = [Vm0, m_0, h_0, n_0, 0, 0];

% Function handles for ODEs
fun_with_NaHT = @(t, Y, I, gap) membranewithNaHT_ODE(t, Y, [I, 1, I, 1, gap, 40]);
fun_without_NaHT = @(t, Y, I, gap) membranewithNaHT_ODE(t, Y, [I, 1, I, 1, gap, 0]);

%% Find the refractory period with NaHT channel
current_gap = initial_gap;
while current_gap <= max_gap && ~AP_triggered_with_NaHT
    [t_with_NaHT, Y_with_NaHT] = ode15s(@(t,Y) fun_with_NaHT(t, Y, I_threshold_unblocked, current_gap), [tstart tend], Y0_with_NaHT);
    Vm_with_NaHT = Y_with_NaHT(:,1);
    if any(Vm_with_NaHT >= V_threshold_with_NaHT)
        refractory_period_with_NaHT = current_gap;
        AP_triggered_with_NaHT = true;
    end
    current_gap = current_gap + gap_increment;
end

%% Find the refractory period without NaHT channel
current_gap = initial_gap;
while current_gap <= max_gap && ~AP_triggered_without_NaHT
    [t_without_NaHT, Y_without_NaHT] = ode15s(@(t,Y) fun_without_NaHT(t, Y, I_threshold_blocked, current_gap), [tstart tend], Y0_without_NaHT);
    Vm_without_NaHT = Y_without_NaHT(:,1);
    if any(Vm_without_NaHT >= V_threshold_without_NaHT)
        refractory_period_without_NaHT = current_gap;
        AP_triggered_without_NaHT = true;
    end
    current_gap = current_gap + gap_increment;
end
% Output the found refractory periods
disp(['Refractory period with NaHT channel: ', num2str(refractory_period_with_NaHT), ' ms']);
disp(['Refractory period without NaHT channel: ', num2str(refractory_period_without_NaHT), ' ms']);


%% Function handles for ODE
fun_blocked = @(t, Y, I, gap) membranewithNaHT_ODE(t, Y, [I, 1, I, 1, gap, 0]);
fun_unblocked = @(t, Y, I, gap) membranewithNaHT_ODE(t, Y, [I, 1, I, 1, gap, 40]);


%% Plot configuration
styles = {'-', '-'}; % Same line style for both gaps
colors = {'k', 'g'}; % Use cell array for colors for more flexibility
lineWidth = 2;       % Set the line width for the plots

%% Blocked NaHT channel plot
figure;
hold on;

% This loop iterates over the specified gap intervals for a neuronal model with the NaHT channel blocked.
for i = 1:length(gap_blocked)
    % `I_threshold_blocked` is the current amplitude set to trigger an action potential, and `gap_blocked(i)` is the current gap interval being tested.
    [t, Y] = ode15s(@(t,Y) fun_blocked(t, Y, I_threshold_blocked, gap_blocked(i)), [tstart tend], Y0_without_NaHT);
    
    % Extract the membrane potential (Vm) from the second column of the solution matrix Y.
    Vm = Y(:,1);
    
    % Plot the membrane potential (Vm) over time (t) with a distinct line style and color for each gap interval.
    plot(t, Vm, 'Color', colors{i}, 'LineStyle', styles{i}, 'LineWidth', lineWidth, ...
         'DisplayName', ['Vm2: t_{gap} = ', num2str(gap_blocked(i)), ' ms']);
end

yline(V_threshold_without_NaHT, 'r--', 'LineWidth', lineWidth, 'DisplayName', 'V_{threshold}');
legend('show', 'Location', 'northeast');
title('Blocked NaHT Channel - Relative Refractory Period');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
hold off;

%% Unblocked NaHT channel plot
figure;
hold on;

% Iterate over gap values for the unblocked NaHT channel
for i = 1:length(gap_unblocked)
    [t, Y] = ode15s(@(t,Y) fun_unblocked(t, Y, I_threshold_unblocked, gap_unblocked(i)), [tstart tend], Y0_with_NaHT);
    Vm = Y(:,1);
    plot(t, Vm, 'Color', colors{i}, 'LineStyle', styles{i}, 'LineWidth', lineWidth, ...
         'DisplayName', ['Vm1: t_{gap} = ', num2str(gap_unblocked(i)), ' ms']);
end

yline(V_threshold_with_NaHT, 'r--', 'LineWidth', lineWidth, 'DisplayName', 'V_{threshold}');
legend('show', 'Location', 'northeast');
title('Unblocked NaHT Channel - Relative Refractory Period');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
hold off;
