% Iterations over different thresholds, tens and cis parameter values.
% Plotting success and max E-level heatmaps shown in Figure 6

% To speed up comment line 71 in
% SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave.mat before 
% running this script to avoid saving simulations movies for each 
% iteration

clear all; close all; clc;

model.ratio = 'rho';
model.zone  = 'Gaussian-like proneural genes';

%% Threshold sweep
thr_values = [0.010 0.013 0.016 0.019];   % detection thresholds
nThr = numel(thr_values);

%% Parameter grids and RNG seeds
cis_exp_values   = linspace(0, 2, 3);
trans_exp_values = linspace(1, 4, 4);
rng_seeds        = 1:30;

nC = numel(cis_exp_values);
nT = numel(trans_exp_values);
nS = numel(rng_seeds);

%% Storage (4D: cis × trans × seed × thr)
success_map  = false(nC, nT, nS, nThr);

numE_map     = nan(nC, nT, nS, nThr);
maxE_NN_map  = nan(nC, nT, nS, nThr);
sdE_map      = nan(nC, nT, nS, nThr);
maxE_SOP_map = nan(nC, nT, nS, nThr);

results = [];

%% Sweep
for iThr = 1:nThr
    thr = thr_values(iThr);

    for iCis = 1:nC
        for iTrans = 1:nT
            for iSeed = 1:nS

                %% --- Model setup ---
                params = SOP_DefaultParams_rho_adapted_weighted(model.zone);
                params.P = 8; params.Q = 8;
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

                params.cis_size_exponent   = cis_exp_values(iCis);
                params.trans_size_exponent = trans_exp_values(iTrans);

                %% --- Geometry randomness ---
                rng(rng_seeds(iSeed));
                params.perimeter = 1 + 0.1.*2.*(rand(k,1) - 0.5);
                params.perimeter(params.SOP.cells) = 0.85;

                %% --- Delamination ---
                params.delam.trigger = 'time';
                params.delam.time = 0;
                params.delam.E_threshold = thr;
                params.delam.SOP_new_perimeter = 0.05;
                params.delam.Ecell_new_perimeter = 1.4;
                params.delam.smooth = true;
                params.delam.duration = 10;

                y0 = SOP_InitialConditions_rho_adapted_weighted(model.zone, params, k);

                %% --- Run simulation ---
                try
                    [yout, tout, params] = ...
                        SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave(model, params, y0);
                catch
                    continue;
                end

                %% Analyze E
                E_max = max(yout(:, 3*k+1:4*k));
                active_E_cells = find(E_max > thr);

                SOP_cells = params.SOP.cells(:);
                NN_cells  = params.SOP.nearest_neighbors(:);
                others    = setdiff(1:k, [SOP_cells; NN_cells]);

                maxE_NN  = max(E_max(NN_cells));
                maxE_SOP = max(E_max(SOP_cells));

                cond1 = ~isempty(intersect(active_E_cells, NN_cells));
                cond2 = all(E_max(SOP_cells) <= thr);
                cond3 = all(E_max(others)    <= thr);
                condSize = isfield(params,'NN_increased') && params.NN_increased;

                condRatio = (maxE_NN / thr) > 10;

                if ~(cond1 && cond2 && cond3 && condSize && condRatio)
                    continue;
                end

                %% --- Metrics ---
                success_map(iCis,iTrans,iSeed,iThr) = true;

                numE      = numel(active_E_cells);
                sdE_cells = std(E_max(active_E_cells));

                numE_map(iCis,iTrans,iSeed,iThr)     = numE;
                maxE_NN_map(iCis,iTrans,iSeed,iThr)  = maxE_NN;
                sdE_map(iCis,iTrans,iSeed,iThr)      = sdE_cells;
                maxE_SOP_map(iCis,iTrans,iSeed,iThr) = maxE_SOP;

                %% Save seed-level
                results = [results; ...
                    cis_exp_values(iCis), ...
                    trans_exp_values(iTrans), ...
                    rng_seeds(iSeed), ...
                    thr, ...
                    numE, maxE_NN, sdE_cells, maxE_SOP];
            end
        end
    end
end

%% Save raw results
results_table = array2table(results, ...
    'VariableNames',{'cis_size_exponent','trans_size_exponent','rng_seed','threshold', ...
                     'num_E_active_cells','maxE_NN', ...
                     'sd_E_between_active_cells','maxE_SOP'});

writetable(results_table,'successful_parameter_metrics_with_thresholds.csv');
disp('Saved successful_parameter_metrics_with_thresholds.csv');

%% Aggregate across RNG seeds AND thresholds
success_prob = mean(success_map, [3 4]);  % average over seed & thr

mean_over_success = @(X) ...
    arrayfun(@(i,j) ...
        mean(X(i,j, success_map(i,j,:,:)), 'omitnan'), ...
        repmat((1:nC)',1,nT), repmat(1:nT,nC,1));

numE_mean     = mean_over_success(numE_map);
maxE_NN_mean  = mean_over_success(maxE_NN_map);
sdE_mean      = mean_over_success(sdE_map);
maxE_SOP_mean = mean_over_success(maxE_SOP_map);

%% Heatmaps
figure1 = figure('Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(cis_exp_values(1:3), trans_exp_values, success_prob(1:3,:)');
axis xy; colorbar; caxis([0 1]);
title('Success probability (avg over thr)');
xlabel('\gamma_{cis}'); ylabel('\gamma_{trans}');
% xlim([0,2])

maps = {numE_mean, maxE_NN_mean, sdE_mean, maxE_SOP_mean};
labels = {'# E-active cells','Max E (NN)','SD E (active cells)','Max E (SOP)'};

for m = 1:4
    nexttile;
    imagesc(cis_exp_values(1:3), trans_exp_values, maps{m}(1:3,:)');
    axis xy; colorbar;
    title([labels{m} ' (avg over thr)']);
    xlabel('\gamma_{cis}'); ylabel('\gamma_{trans}');
    % xticks([0,1,2])
    yticks([1,2,3,4])
    title(labels{m});
    colormap('sky')
    if m==2
       clim([0,Inf])
    end
    % xlim([0,2])
end

sgtitle('cis × trans heatmaps averaged over RNG & threshold');
saveas(figure1,'RobustnessAveragedHeatmaps_withThresholds.pdf');
saveas(figure1,'RobustnessAveragedHeatmaps_withThresholds.png');