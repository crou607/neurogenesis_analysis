% Iterations over different delamination durations

% To speed up comment line 71 in
% SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave.mat before 
% running this script to avoid saving simulations movies for each 
% iteration


clear all; close all; clc;

model.ratio = 'rho';
model.zone  = 'Gaussian-like proneural genes';

thr = 0.012;

%% Delamination durations and RNG seeds
delam_durations = [5 10 15];
rng_seeds       = 1:100;

nD = numel(delam_durations);
nS = numel(rng_seeds);

%% Storage
success_map   = false(nD, nS);
timeSS_map    = nan(nD, nS);      % time to steady state
maxE_NN_map   = nan(nD, nS);      % maximum E in NN cells

results = [];

%% Sweep
for iD = 1:nD
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

        %% Fixed exponents
        params.cis_size_exponent   = 1;
        params.trans_size_exponent = 3;

        %% Geometry randomness
        rng(rng_seeds(iSeed));
        params.perimeter = 1 + 0.1.*2.*(rand(k,1) - 0.5);
        params.perimeter(params.SOP.cells) = 0.85;

        %% Delamination
        params.delam.trigger = 'time';
        params.delam.time = 0;
        params.delam.E_threshold = thr;
        params.delam.SOP_new_perimeter = 0.05;
        params.delam.Ecell_new_perimeter = 1.4;
        params.delam.smooth = true;
        params.delam.duration = delam_durations(iD);

        y0 = SOP_InitialConditions_rho_adapted_weighted(model.zone, params, k);

        %% Run simulation
        try
            [yout, tout, params] = ...
                SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave(model, params, y0);
        catch
            continue;
        end

        %% Analyze E dynamics
        E_all = yout(:, 3*k+1:4*k);        % E for all cells over time
        E_max = max(E_all);                % peak E per cell

        SOP_cells = params.SOP.cells(:);
        NN_cells  = params.SOP.nearest_neighbors(:);
        others    = setdiff(1:k, [SOP_cells; NN_cells]);

        active_E_cells = find(E_max > thr);

        maxE_NN = max(E_max(NN_cells));

        %% Success conditions 
        cond1 = ~isempty(intersect(active_E_cells, NN_cells));
        cond2 = all(E_max(SOP_cells) <= thr);
        cond3 = all(E_max(others)    <= thr);
        condSize  = isfield(params,'NN_increased') && params.NN_increased;
        condRatio = (maxE_NN / thr) > 10;

        if ~(cond1 && cond2 && cond3 && condSize && condRatio)
            continue;
        end

        %% Compute time to steady state (robust plateau detection)
        mean_NN_E = mean(E_all(:, NN_cells), 2);
        final_value = mean_NN_E(end);

        tol = 0.02 * final_value;   % 2% tolerance band
        within_band = abs(mean_NN_E - final_value) < tol;

        steady_idx = NaN;

        for i = 1:length(mean_NN_E)
            if all(within_band(i:end))
                steady_idx = i;
                break
            end
        end

        if isnan(steady_idx)
            continue;   % never reached steady state
        end

        time_to_ss = tout(steady_idx);

        %% Store results
        success_map(iD,iSeed) = true;
        timeSS_map(iD,iSeed)  = time_to_ss;
        maxE_NN_map(iD,iSeed) = maxE_NN;

        results = [results; ...
            delam_durations(iD), ...
            rng_seeds(iSeed), ...
            time_to_ss, ...
            maxE_NN];
    end
end

%% Save CSV
results_table = array2table(results, ...
    'VariableNames',{'delamination_duration', ...
                     'rng_seed', ...
                     'time_to_steady_state', ...
                     'maxE_NN'});

writetable(results_table,'time_to_steady_state_and_maxE_successful.csv');
disp('Saved time_to_steady_state_and_maxE_successful.csv');

