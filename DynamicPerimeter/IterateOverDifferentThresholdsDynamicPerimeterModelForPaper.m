% Iterations over different thresholds, kappa_p and NB initial perimeters.
% Plotting success and max E-level heatmaps shown in Figure S4

% To speed up comment line 79 in
% SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave.mat before 
% running this script to avoid saving simulations movies for each 
% iteration


%% 3D SWEEP: perimeter_kappa × SOP p0 × threshold
clear; close all; clc;

model.ratio = 'rho';
model.zone  = 'Gaussian-like proneural genes';

%% Parameter grids
perimeter_kappa_vals = [0.15 0.2 0.25];
SOP_p0_vals          = linspace(0.65, 0.85, 5);
thr_values           = [0.015 0.02 0.025];   % can expand again later
rng_seeds            = 1:10;

%% Storage
Results = [];
row = 0;

%% Sweep
for thr = thr_values
for k_del = perimeter_kappa_vals
for p0 = SOP_p0_vals
for seed = rng_seeds

    row = row + 1;

    %% Model setup
    params = SOP_DefaultParams_rho_adapted_weighted(model.zone);
    params.P = 8; 
    params.Q = 8;
    params.Tmax = 10;
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
    params.cis_size_exponent   = 0;
    params.trans_size_exponent = 3;

    %% Perimeter dynamics 
    params.perimeter_kappa       = k_del;
    params.perimeter_hill        = 2;
    params.perimeter_E_hill      = 1;
    params.perimeter_theta_delta = 0.4;
    params.perimeter_theta_E     = 0.4;

    %% Geometry
    rng(seed);
    params.perimeter = 1 + 0.1.*2.*(rand(k,1) - 0.5);

    SOP_cells = params.SOP.cells(:);
    NN_cells  = params.SOP.nearest_neighbors(:);
    others    = setdiff(1:k, [SOP_cells; NN_cells]);

    params.perimeter(SOP_cells) = p0;

    %% Initial conditions 
    y0 = SOP_InitialConditions_rho_adapted_weighted(model.zone, params, k);

    %% Run simulation 
    success = false;
    numE = NaN;
    maxE_NN = NaN;
    maxE_SOP = NaN;
    finalNB = NaN;
    finalNN = NaN;
    finalOther = NaN;

    try
        [yout, tout, params] = ...
            SOP_multicell_LI_adapted_DYNAMIC_perimeter_HillsNoRecovery(model, params, y0);
    catch
        Results(row,:) = [k_del p0 thr seed false ...
                          numE maxE_NN maxE_SOP ...
                          finalNB finalNN finalOther];
        continue;
    end

    %% Analyze E
    E_traj = yout(:, 3*k+1:4*k);
    E_max  = max(E_traj);

    active_E_cells = find(E_max > thr);

    maxE_NN  = max(E_max(NN_cells));
    maxE_SOP = max(E_max(SOP_cells));

    %% Analyze perimeter
    p_traj  = yout(:, 5*k+1:6*k);
    p_final = p_traj(end,:);

    NB_min_perimeter = min(min(p_traj(:, SOP_cells)));

    finalNB    = mean(p_final(SOP_cells));
    finalNN    = mean(p_final(NN_cells));
    finalOther = mean(p_final(others));

    %% Success conditions
    cond1 = ~isempty(intersect(active_E_cells, NN_cells));
    cond2 = all(E_max(SOP_cells) <= thr);
    cond3 = all(E_max(others)    <= thr);
    condRatio = (maxE_NN / thr) > 10;
    condSizeFlag = isfield(params,'NN_increased') && params.NN_increased;

    condNBshrink   = NB_min_perimeter <= 0.3;
    condOthersStay = finalOther >= 0.6;

    success = all([cond1 cond2 cond3 condRatio ...
                   condSizeFlag condNBshrink condOthersStay]);

    %% Store
    numE = numel(active_E_cells);

    Results(row,:) = [ ...
        k_del, p0, thr, seed, ...
        success, numE, ...
        maxE_NN, maxE_SOP, ...
        finalNB, finalNN, finalOther];

end
end
end
end

%% ===============================================================
% Convert to table
%% ===============================================================
Results = array2table(Results, ...
    'VariableNames',{ ...
    'kappa_perimeter', ...
    'SOP_initial_perimeter', ...
    'threshold', ...
    'rng_seed', ...
    'success', ...
    'num_E_active', ...
    'maxE_NN', ...
    'maxE_SOP', ...
    'final_NB_perimeter', ...
    'final_NN_perimeter', ...
    'final_other_perimeter'});

%% HEATMAPS
nK = length(perimeter_kappa_vals);
nP = length(SOP_p0_vals);

success_map  = nan(nK,nP);

% These will use only successful simulations
maxE_NN_map  = nan(nK,nP);
maxE_SOP_map = nan(nK,nP);
NB_map       = nan(nK,nP);
NN_map       = nan(nK,nP);
Other_map    = nan(nK,nP);

SuccessOnly = Results(Results.success == 1, :);

for i = 1:nK
    for j = 1:nP

        %% SUCCESS PROBABILITY → ALL simulations
        mask_all = Results.kappa_perimeter == perimeter_kappa_vals(i) & ...
                   Results.SOP_initial_perimeter == SOP_p0_vals(j);

        success_map(i,j) = mean(Results.success(mask_all));

        %% OTHER METRICS → SUCCESSFUL simulations only
        mask_succ = SuccessOnly.kappa_perimeter == perimeter_kappa_vals(i) & ...
                    SuccessOnly.SOP_initial_perimeter == SOP_p0_vals(j);

        if any(mask_succ)

            maxE_NN_map(i,j)  = mean(SuccessOnly.maxE_NN(mask_succ));
            maxE_SOP_map(i,j) = mean(SuccessOnly.maxE_SOP(mask_succ));
            NB_map(i,j)       = mean(SuccessOnly.final_NB_perimeter(mask_succ));
            NN_map(i,j)       = mean(SuccessOnly.final_NN_perimeter(mask_succ));
            Other_map(i,j)    = mean(SuccessOnly.final_other_perimeter(mask_succ));

        end
    end
end
