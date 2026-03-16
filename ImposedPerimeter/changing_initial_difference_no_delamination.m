% Iterations over various initial NB-other cell size differences without 
% delamination. Plotting success probability for each size difference 
% shown in Figure 6.

% To speed up comment line 71 in
% SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave.mat before 
% running this script to avoid saving simulations movies for each 
% iteration

clear all; close all; clc;

model.ratio = 'rho';
model.zone  = 'Gaussian-like proneural genes';

thr = 0.012;

%% Sweep over initial SOP perimeter (p0)
SOP_p0_values = linspace(0.1, 0.9, 17);   % adjust range if needed
rng_seeds     = 1:30;

nP = numel(SOP_p0_values);
nS = numel(rng_seeds);

%% Storage
success_map    = false(nP, nS);

maxE_NB_map    = nan(nP, nS);   % max E in SOP (NB)
maxE_NN_map    = nan(nP, nS);   % max E in NN
sdE_NN_map     = nan(nP, nS);   % SD of E across NN cells (at peak)

results = [];

%% Sweep
for iP = 1:nP
    for iSeed = 1:nS

        %% Model setup 
        params = SOP_DefaultParams_rho_adapted_weighted(model.zone);
        params.P = 8;
        params.Q = 8;
        params.Tmax = 20;
        k = params.P * params.Q;

        params.beta.d = 1.62;
        params.beta.E = 1.62;
        params.beta.n = 1.52;
        params.beta.A = 1.62;
        params.T.s = 0.8;
        params.T.E = 0.0051;
        params.c.s = 3;
        params.c.E = 3;
        params.Kappa.t = 0.8;
        params.Kappa.InhCombined = 0.4;

        params.cis_size_exponent   = 1;
        params.trans_size_exponent = 3;

        %% Geometry randomness
        rng(rng_seeds(iSeed));
        params.perimeter = 1 + 0.1.*2.*(rand(k,1) - 0.5);

        SOP_cells = params.SOP.cells(:);
        NN_cells  = params.SOP.nearest_neighbors(:);

        %% Set SOP initial perimeter (p0 sweep)
        params.perimeter(SOP_cells) = SOP_p0_values(iP);

        %% --- Disable delamination completely ---
        %% --- Delamination ---
        params.delam.trigger = 'time';
        params.delam.time = 0;
        params.delam.E_threshold = thr;
        params.delam.SOP_new_perimeter = SOP_p0_values(iP);
        params.delam.Ecell_new_perimeter = 1.4;
        params.delam.smooth = true;
        params.delam.duration = 10;
        

        y0 = SOP_InitialConditions_rho_adapted_weighted(model.zone, params, k);

        %% --- Run simulation ---
        try
            [yout, ~, params] = ...
                SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave(model, params, y0);
        catch
            continue;
        end

        %% Analyze E dynamics
        E_all = yout(:, 3*k+1:4*k);    % E over time
        E_max = max(E_all);           % peak E per cell

        others = setdiff(1:k, [SOP_cells; NN_cells]);

        active_E_cells = find(E_max > thr);

        maxE_NB = max(E_max(SOP_cells));
        maxE_NN = max(E_max(NN_cells));
        sdE_NN  = std(E_max(NN_cells));

        %% Success conditions
        cond1 = ~isempty(intersect(active_E_cells, NN_cells));
        cond2 = all(E_max(SOP_cells) <= thr);
        cond3 = all(E_max(others)    <= thr);
        condSize  = isfield(params,'NN_increased') && params.NN_increased;
        condRatio = (maxE_NN / thr) > 10;

        if ~(cond1 && cond2 && cond3 && condSize && condRatio)
            continue;
        end

        %% Store successful results
        success_map(iP,iSeed)  = true;

        maxE_NB_map(iP,iSeed)  = maxE_NB;
        maxE_NN_map(iP,iSeed)  = maxE_NN;
        sdE_NN_map(iP,iSeed)   = sdE_NN;

        results = [results; ...
            SOP_p0_values(iP), ...
            rng_seeds(iSeed), ...
            maxE_NB, ...
            maxE_NN, ...
            sdE_NN];
    end
end

%% Save CSV (successful runs only)
results_table = array2table(results, ...
    'VariableNames',{'SOP_initial_perimeter', ...
                     'rng_seed', ...
                     'maxE_NB', ...
                     'maxE_NN', ...
                     'sdE_NN'});

writetable(results_table,'SOP_p0_sweep_successful.csv');
disp('Saved SOP_p0_sweep_successful.csv');

%% Success probability per p0
success_prob = mean(success_map,2);

disp('Success probability per SOP p0:')
disp(table(SOP_p0_values', success_prob, ...
    'VariableNames',{'SOP_initial_perimeter','success_probability'}))

%% Aggregate statistics over successful runs
mean_over_success = @(X) ...
    arrayfun(@(i) ...
        mean(X(i, success_map(i,:)), 'omitnan'), ...
        1:nP)';

std_over_success = @(X) ...
    arrayfun(@(i) ...
        std(X(i, success_map(i,:)), 'omitnan'), ...
        1:nP)';

maxE_NB_mean = mean_over_success(maxE_NB_map);
maxE_NN_mean = mean_over_success(maxE_NN_map);
sdE_NN_mean  = mean_over_success(sdE_NN_map);

maxE_NB_std = std_over_success(maxE_NB_map);
maxE_NN_std = std_over_success(maxE_NN_map);
sdE_NN_std  = std_over_success(sdE_NN_map);

%% Error plots
fig1 = figure('Color','w');
plot(SOP_p0_values(2:end), success_prob(2:end), 'LineWidth',2);
xlim([0.15 0.9])
ylim([0 1]);
xlabel('SOP initial perimeter (p0)');
ylabel('Success probability');
title('Success probability');
saveas(fig1,'probability_p0_plot.pdf');
saveas(fig1,'probability_p0_plot.fig');


fig2 = figure('Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
% --- Max E in NB ---
nexttile;
errorbar(SOP_p0_values, maxE_NB_mean, maxE_NB_std, ...
    'LineWidth',1.5);
xlabel('SOP initial perimeter (p0)');
ylabel('Max E (NB)');
title('NB peak E');

% --- Max E in NN ---
nexttile;
errorbar(SOP_p0_values, maxE_NN_mean, maxE_NN_std, ...
    'LineWidth',1.5);
xlabel('SOP initial perimeter (p0)');
ylabel('Max E (NN)');
title('NN peak E');

% --- SD E in NN ---
nexttile;
errorbar(SOP_p0_values, sdE_NN_mean, sdE_NN_std, ...
    'LineWidth',1.5);
xlabel('SOP initial perimeter (p0)');
ylabel('SD E (NN)');
title('NN variability');

sgtitle('SOP p0 sweep – mean ± SD across successful runs');

saveas(fig2,'SOP_p0_errorplots.pdf');
saveas(fig2,'SOP_p0_errorplots.fig');