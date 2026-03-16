% Script that runs imposed perimeter model shown in Figure 6 for a specific
% parameter choice

% Adapted from Troost et al., 2023

%% RUN SCRIPT for SOP Weighted Lateral Inhibition Model

clear; close all; clc;

%% Model setup

model.ratio = 'rho';
model.zone  = 'Gaussian-like proneural genes';

%% Simulation parameters
params = SOP_DefaultParams_rho_adapted_weighted(model.zone);

% Lattice geometry
params.P = 8;   % rows
params.Q = 8;   % columns
params.Tmax = 20;  % simulation duration

% Control SOP patterning 
% params.nSOP = 2;        % number of SOPs
% params.minDist = 3;     % minimum distance (in graph hops)

% Parameter scaling 
params.beta.d = 1.62;   % base Delta production
params.beta.E = 1.62;   % base E(spl) production     %%% NOTE %%%
params.beta.n = 1.52;   % base Notch production

% Transcription and thresholds 
params.T.E = 0.0051; % Hill function threshold for inhibitng activator expression by E(spl)  
% params.c.s = 2; % Hill coefficient for promoting E(spl) expression by trans-activation signal  %%% NOTE %%%
params.c.E = 3; % Hill coefficient for inhibitng activator expression by E(spl)                %%% NOTE %%%
params.Kappa.t = 0.8;            % transactivation rate
params.Kappa.InhCombined = 0.4;  % base cis interaction scaling                                %%% NOTE %%%

k = params.P * params.Q;

% Perimeter and delamination/enhancement       
rng(21)
params.perimeter = 1 + 0.1.*2.*(rand(k,1) - 0.5);
params.perimeter(params.SOP.cells) = 0.85;  % smaller center cell


params.delam.trigger = 'time';               % (default: no event)
params.delam.time = 0;     %                (model time units)  -- used if trigger='time'
params.delam.E_threshold = 0.015;  %           (used if trigger='E', compares E > threshold)
params.delam.SOP_new_perimeter = 0.05;
params.delam.Ecell_new_perimeter = 1.4;  %% Note: this is fold change (multiplied with p0) so the differences maintain
params.delam.smooth = true; %true;
params.delam.duration = 10;  

params.beta.A = 1.62;  % base Activator production rate (1.62 in paper)
params.c.s = 3; % Hill coefficient for promoting E(spl) expression by trans-activation signal
params.T.s = 0.8; % Hill function threshold for promoting E(spl) expression by trans-activation signal
params.cis_size_exponent = 1;  % how much smaller cells increase cis inhibition
params.trans_size_exponent = 3; % trans interaction scaling

%% INITIAL CONDITIONS 
y0 = SOP_InitialConditions_rho_adapted_weighted(model.zone, params, k);
% 
%% RUN SIMULATION 
fprintf('\nRunning SOP weighted lateral inhibition simulation...\n\n');
[yout, tout, params] = SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave(model, params, y0);

%% Map model time to real seconds
real_time_per_model_unit = 150;  % 1 model unit = 150s

t_real = tout * real_time_per_model_unit;

% Optional: resample every 30 s (like your experimental frames)
frame_interval = 150;  % seconds
t_uniform = 0:frame_interval:max(t_real);
y_uniform = interp1(t_real, yout, t_uniform, 'pchip');


%%
figure(10); clf; hold on;

tmin = t_uniform/60;
SOPcells = params.SOP.cells(:);

% Extract blocks
D_SOP   = y_uniform(:, SOPcells);                % Delta
DA_SOP  = y_uniform(:, k + SOPcells);            % Activated Delta
N_SOP   = y_uniform(:, 2*k + SOPcells);          % Notch
E_SOP   = y_uniform(:, 3*k + SOPcells);          % E(spl)
Aac_SOP = y_uniform(:, 4*k + SOPcells);          % Activator

% Colors
cD  = [0 0.4 1];
cDA = [0.1 0.7 1];
cN  = [0.8 0.2 0];
cE  = [0.6 0 0.8];
cA  = [0 0.6 0];

% Plot mean ± SEM
shadedMeanPlot(tmin, D_SOP,  cD,  'Delta');
shadedMeanPlot(tmin, DA_SOP, cDA, 'Activated Delta');
shadedMeanPlot(tmin, N_SOP,  cN,  'Notch');
shadedMeanPlot(tmin, E_SOP,  cE,  'E(spl)');
shadedMeanPlot(tmin, Aac_SOP, cA,  'Activator');

legend('Location','best');
xlabel('Time (min)');
ylabel('Level (mean ± SD)');
title('SOP cells');
xlim([0 max(tmin)]);
ylim([0 1.8]);
saveas(figure(10),'NB_levels.pdf')
%%
figure(12); clf; hold on;
NNcells = params.SOP.nearest_neighbors(:);


% Extract blocks
D_NN   = y_uniform(:, NNcells);
DA_NN  = y_uniform(:, k + NNcells);
N_NN   = y_uniform(:, 2*k + NNcells);
E_NN   = y_uniform(:, 3*k + NNcells);
Aac_NN = y_uniform(:, 4*k + NNcells);

% Plot mean ± SEM
shadedMeanPlot(tmin, D_NN,  cD,  'Delta');
shadedMeanPlot(tmin, DA_NN, cDA, 'Activated Delta');
shadedMeanPlot(tmin, N_NN,  cN,  'Notch');
shadedMeanPlot(tmin, E_NN,  cE,  'E(spl)');
shadedMeanPlot(tmin, Aac_NN, cA,  'Activator');

legend('Location','best');
xlabel('Time (min)');
ylabel('Level (mean ± SD)');
title('Nearest neighbour cells');
xlim([0 max(tmin)]);
ylim([0 1.8]);
saveas(figure(12),'NN_levels.pdf')
%%
figure(14); clf; hold on;

% Plot mean ± SEM
shadedMeanPlot(tmin, E_SOP,  'g',  'E(spl) SOP');
shadedMeanPlot(tmin, E_NN, 'b',  'E(spl) NN');

legend('Location','best');
xlabel('Time (min)');
ylabel('Level (mean ± SD)');
title('E(spl) levels');
xlim([0 max(tmin)]);
ylim([0 0.15]);
% hold on
% x1 =[5 5];
% y1 = [0 0.25];
% plot(x1,y1,'--k')
saveas(figure(14),'E_levels.pdf')

figure; hold on
yyaxis right
shadedMeanPlot(tmin, E_SOP,  'g',  'E(spl) SOP');
shadedMeanPlot(tmin, E_NN, 'b',  'E(spl) NN');
ylim([0 0.18])

yyaxis left
ylim([0 1.6])
saveas(figure,'E_levels_and_perimeters.pdf')
