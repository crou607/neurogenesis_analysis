% Script with default model parameters

% Adapted from Troost et al., 2023

function params = SOP_DefaultParams_rho_adapted_weighted(model)
%
% - NOTE: Perimeter-based weighting and many behaviors are implemented in
%   SOP_multicell_LI_adapted_weighted (runtime), but this file prepares
%   SOP indices and default parameter values.

params.Tmax = 100;
params.P = 8;
params.Q = 8;

% Base production rates (dimensionless)
params.beta.d = 1*1.62;  % base delta production
params.beta.n = 1*1.52;  % base notch production
params.beta.E = 1*1.62;  % base E(spl) production
params.beta.A = 1*1.62;  % NEW: base activator production (used in adapted model)

% Activation/inactivation parameters
params.alpha.plus = 0.6; % NEW: activation rate of d_A by A (role similar to alpha_N previously)
params.alpha.minus = 0.4; % retention/decay related parameter

% Cis / trans inhibition parameters
params.Kappa.t = 1;               % trans-activation rate
params.Kappa.InhCombined = 1;     % NEW: combined cis-inhibition scaling (base value)
params.cis_size_exponent = 1;     % NEW: exponent controlling sensitivity of cis to cell size (perimeter)
params.trans_size_exponent = 1;

% Thresholds & Hill coefficients (for transcription / nonlinear response)
params.T.E = 0.0051; % E threshold used for A production repression (used in Hill function)
params.T.s = 5;      % signal threshold for E transcription (used in Hill function)
params.c.s = 2;      % Hill coefficient for signal→E
params.c.E = 3;      % Hill coefficient for E→A (or similar) nonlinearities

% (old parameters removed or renamed in this adapted version)

switch model
    case 'Mib1 mutual inhibition zone' 
        % nothing special here in adapted version (kept for API stability)
        
    case 'Neur lateral inhibition zone'
        % nothing special here in adapted version (kept for API stability)
        
    case 'Gaussian-like proneural genes'
        % Create multiple SOPs (Gaussian-like proneural gene pattern)
 
        Gaussian_steps = [1 0.61 0.14]; % normalized amplitudes for center, nearest, next-nearest
                                         % (typical "Gaussian" discretization over hex lattice)
                
        C = getconnectivityM(params.P, params.Q); % adjacency (unweighted) matrix
        k = params.P * params.Q;                   % number of cells in lattice
        
        % SOP selection: FIXED in this file to make runs reproducible
        % Previously we had a stochastic pick with distance constraints.
        % Here we set SOP_cells to fixed indices so multiple runs pick the same SOPs.
        % (You can change the indices [21,42] to any cells that are far apart for your P/Q.)
        SOP_cells = [21, 42, 46];
        
        % build neighborhoods for each SOP (nearest & next-nearest)
        SOP_nearest_neighbors = [];
        SOP_next_nearest_neighbors = [];
        for s = SOP_cells
            % Note: the adjacency matrix built here uses w=[1/6 ...] so neighbors show 1/6.
            % We therefore find entries equal to 1/6 (the neighbor weights).
            [~, nn]  = find(C(s,:) == 1/6);      % nearest neighbors of SOP
            [~, nnn] = find(C(nn',:) == 1/6);    % neighbors of those neighbors
            nnn = setdiff(nnn, [s ; nn']);      % remove the SOP and nearest neighbors -> next-nearest
            SOP_nearest_neighbors = [SOP_nearest_neighbors, nn];
            SOP_next_nearest_neighbors = [SOP_next_nearest_neighbors, nnn];
        end
        
        % Deduplicate indices (in case neighborhoods overlap)
        SOP_nearest_neighbors = unique(SOP_nearest_neighbors);
        SOP_next_nearest_neighbors = setdiff(unique(SOP_next_nearest_neighbors), ...
                                             [SOP_cells, SOP_nearest_neighbors]);
        
        % Save chosen SOP indices in params so other functions (initial conditions,
        % visualization) can use the same SOP set.
        params.SOP.cells = SOP_cells;
        params.SOP.nearest_neighbors = SOP_nearest_neighbors;
        params.SOP.next_nearest_neighbors = SOP_next_nearest_neighbors;

        
    otherwise
        warning('Unexpected model type. See SOP_defaultparams description')
end


end


%% additional functions (unchanged connectivity helpers)
function C=getconnectivityM(P,Q)

k=P*Q; %number of cells
C=zeros(k,k); %This is the connectivity matrix
w=[1/6 1/6 1/6 1/6 1/6 1/6]; % uniform weight for the 6 neighbors
 
% calculating the connectivity matrix for hex lattice
for s=1:k
    kneighbour=findneighbourhex(s,P,Q); %finds the 6 neighbors of cell s
    for r=1:6
        C(s,kneighbour(r))=w(r);
    end
end
end

function out = findneighbourhex(ind,P,Q)
% Find 6 neighbors (hex grid) using column offset coordinates
[p,q] = ind2pq(ind,P);

% above and below:
out(1) = pq2ind(mod(p,P)+1,q,P);
out(2) = pq2ind(mod(p-2,P)+1,q,P);

% left & right column indices (with wrap-around)
qleft = mod(q-2,Q)+1;
qright = mod(q,Q)+1;

% hex offset depends on column parity:
if mod(q,2) == 1   % odd column
    pup   = p;
    pdown = mod(p-2,P)+1;
else               % even column
    pup   = mod(p, P)+1;
    pdown = p;
end

% neighbors in left & right columns (two each)
out(3) = pq2ind(pup,qleft,P);
out(4) = pq2ind(pdown,qleft,P);
out(5) = pq2ind(pup,qright,P);
out(6) = pq2ind(pdown,qright,P);
end

function ind=pq2ind(p,q, P)
% Convert (p,q) lattice coordinates to linear index
ind = p + (q-1)*P;
end

function [p,q]=ind2pq(ind, P)
% Convert linear index to (p,q)
q = 1+floor((ind-1)/P);
p = ind - (q-1)*P;
end