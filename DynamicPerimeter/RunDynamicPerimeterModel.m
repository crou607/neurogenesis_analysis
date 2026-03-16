% Script that runs dynamic perimeter model shown in Figure S4 for a 
% specific parameter choice

% Adapted from Troost et al., 2023

%% RUN SCRIPT for SOP Weighted Lateral Inhibition Model
%% Dynamic perimeter driven by Activated Delta

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

% Parameter scaling 
params.beta.d = 1.62;
params.beta.E = 1.62;
params.beta.n = 1.52;

% Transcription and thresholds 
params.T.E = 0.0051;
params.c.E = 3;
params.Kappa.t = 0.8;
params.Kappa.InhCombined = 0.4;

params.beta.A = 1.62;
params.c.s = 3;
params.T.s = 0.8;

params.cis_size_exponent   = 0;
params.trans_size_exponent = 3;

k = params.P * params.Q;

%% ===== Initial perimeter (baseline geometry) =====
% rng(17)
rng(8)
params.perimeter = 1 + 0.1.*2.*(rand(k,1) - 0.5);
params.perimeter(params.SOP.cells) = 0.75;   % smaller SOPs initially

%% ===== NEW: Perimeter dynamics parameters =====
params.perimeter_kappa = 0.2;    % shrink by d_A
params.perimeter_hill = 2;       % Hill exponent for d_A
params.perimeter_E_hill = 1;      % Hill exponent for E

params.perimeter_theta_delta = 0.4; % Hill threshold for d_A
params.perimeter_theta_E = 0.4;     % Hill threshold for E

% params.perimeter_t_on = 5;   % e.g. 10
% params.perimeter_t_slope  = 10;
%% INITIAL CONDITIONS (unchanged)
y0 = SOP_InitialConditions_rho_adapted_weighted(model.zone, params, k);

%% RUN SIMULATION
fprintf('\nRunning SOP weighted lateral inhibition simulation (dynamic perimeter)...\n\n');
[yout, tout, params] = SOP_multicell_LI_adapted_DYNAMIC_perimeter_HillsNoRecovery(model, params, y0);

%% Map model time to real seconds
real_time_per_model_unit = 150;  % 1 model unit = 150 s
t_real = tout * real_time_per_model_unit;

% Resample every 30 s
frame_interval = 150;
t_uniform = 0:frame_interval:max(t_real);
y_uniform = interp1(t_real, yout, t_uniform, 'pchip');

tmin = t_uniform / 60;

%% ===== STATE BLOCKS =====
D   = y_uniform(:, 1:k);
DA  = y_uniform(:, k+1:2*k);
N   = y_uniform(:, 2*k+1:3*k);
E   = y_uniform(:, 3*k+1:4*k);
Aac = y_uniform(:, 4*k+1:5*k);
P   = y_uniform(:, 5*k+1:6*k);   % <-- perimeter (NEW)

SOPcells = params.SOP.cells(:);
NNcells  = params.SOP.nearest_neighbors(:);

%% ===== COLORS =====
cD  = [0 0.4 1];
cDA = [0.1 0.7 1];
cN  = [0.8 0.2 0];
cE  = [0.6 0 0.8];
cA  = [0 0.6 0];
cP  = [0.3 0.3 0.3];

%% ===== SOP CELLS =====
figure(10); clf; hold on;

shadedMeanPlot(tmin, D(:,SOPcells),   cD,  'Delta');
shadedMeanPlot(tmin, DA(:,SOPcells),  cDA, 'Activated Delta');
shadedMeanPlot(tmin, N(:,SOPcells),   cN,  'Notch');
shadedMeanPlot(tmin, E(:,SOPcells),   cE,  'E(spl)');
shadedMeanPlot(tmin, Aac(:,SOPcells), cA,  'Activator');

legend('Location','best');
xlabel('Time (min)');
ylabel('Level (mean ± SD)');
title('SOP cells');
xlim([0 max(tmin)]);
ylim([0 1.8]);
saveas(gcf,'NB_levels.pdf')

%% ===== NEAREST NEIGHBORS =====
figure(12); clf; hold on;

shadedMeanPlot(tmin, D(:,NNcells),   cD,  'Delta');
shadedMeanPlot(tmin, DA(:,NNcells),  cDA, 'Activated Delta');
shadedMeanPlot(tmin, N(:,NNcells),   cN,  'Notch');
shadedMeanPlot(tmin, E(:,NNcells),   cE,  'E(spl)');
shadedMeanPlot(tmin, Aac(:,NNcells), cA,  'Activator');

legend('Location','best');
xlabel('Time (min)');
ylabel('Level (mean ± SD)');
title('Nearest neighbour cells');
xlim([0 max(tmin)]);
ylim([0 1.8]);
saveas(gcf,'NN_levels.pdf')

%% ===== E(spl) comparison =====
figure(14); clf; hold on;

shadedMeanPlot(tmin, E(:,SOPcells), 'g', 'E(spl) SOP');
shadedMeanPlot(tmin, E(:,NNcells),  'b', 'E(spl) NN');

legend('Location','best');
xlabel('Time (min)');
ylabel('Level (mean ± SD)');
title('E(spl) levels');
xlim([0 max(tmin)]);
ylim([0 0.6]);
saveas(gcf,'E_levels.pdf')

%% ===== NEW: PERIMETER DYNAMICS =====
figure(16); clf; hold on;

shadedMeanPlot(tmin, P(:,SOPcells), cP, 'Perimeter SOP');
shadedMeanPlot(tmin, P(:,NNcells),  [0.6 0.2 0.2], 'Perimeter NN');

legend('Location','best');
xlabel('Time (min)');
ylabel('Perimeter (a.u.)');
title('Emergent perimeter dynamics');
xlim([0 max(tmin)]);
saveas(gcf,'Perimeter_levels.pdf')



fig0 = figure; hold on
yyaxis right
shadedMeanPlot(tmin, E(:,SOPcells), 'g', 'E(spl) SOP');
shadedMeanPlot(tmin, E(:,NNcells),  'b', 'E(spl) NN');
ylim([0 0.7])

yyaxis left
ylim([0 1.6])
shadedMeanPlot(tmin, P(:,SOPcells), cP, 'Perimeter SOP');
shadedMeanPlot(tmin, P(:,NNcells),  [0.6 0.2 0.2], 'Perimeter NN');

saveas(fig0,'E_levels_and_perimeters.fig')







