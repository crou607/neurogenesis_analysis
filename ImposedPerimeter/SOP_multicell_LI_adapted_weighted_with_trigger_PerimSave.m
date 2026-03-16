% Main script with ODEs and solvers for imposed perimeter model

% Adapted from Troost et al., 2023

function [yout,tout,params] = SOP_multicell_LI_adapted_weighted_with_trigger_PerimSave(model,params,y0)

%% prepare
pp = genpath('functions');
addpath(pp)

% load defaults if absent
if(nargin < 2)
    switch model.ratio
        case 'rho'
            params = SOP_DefaultParams_rho_adapted(model.zone);
        case 'epsilon'
            params = SOP_DefaultParams_epsilon_adapted(model.zone);
    end
end

Tmax = params.Tmax; 
tspan=[0 Tmax];

% Base connectivity matrix (unweighted neighbor matrix)
C_base = getconnectivityM(params.P,params.Q);
params.connectivity = C_base;
k = params.P * params.Q;

% Ensure perimeter baseline exists
if ~isfield(params,'perimeter') || numel(params.perimeter)~=k
    params.perimeter = ones(k,1);  % default uniform
end
params.perimeter0 = params.perimeter(:);  % baseline (initial) perimeters

% Precompute static things used by li
params.k = k;
params.C_base = C_base;

% If y0 not provided, compute
if(nargin < 3)
    switch model.ratio
        case 'rho'
            y0 = SOP_InitialConditions_rho_adapted_weighted(model.zone,params,k);
        case 'epsilon'
            y0 = SOP_InitialConditions_epsilon_adapted(model.zone,params,k);
    end
end

clear_global_perimeter_store();

%% Run solver (li uses globals to store perimeter history)
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[tout,yout] = ode15s(@li,tspan,y0,opts,params);

%% Copy global store into params and clear globals
global PERIM_STORE T_STORE P_COUNTER

if exist('PERIM_STORE','var') && ~isempty(PERIM_STORE)
    % PERIM_STORE is k x N
    nrec = P_COUNTER - 1;
    params.perimeter_store = PERIM_STORE(:,1:nrec);
    params.perimeter_time  = T_STORE(1:nrec);
else
    params.perimeter_store = params.perimeter0(:);
    params.perimeter_time = tout(:)';
end
% clear globals
clear_global_perimeter_store();

%% Make movies and plots
F  = movielattice_and_perimeter_movies(tout,yout,k,params);

%% POST-PROCESS: Detect whether any NN cell increased in perimeter
% default output (in case no delamination or no data)
params.NN_increased = false;
params.dp = [];
% need perimeter history
if isfield(params,'perimeter_store') && isfield(params,'perimeter_time')

    perim_store = params.perimeter_store;   % k x N

    params.dp = perim_store(:,end) - perim_store(:,1);

    NNcells = params.SOP.nearest_neighbors(:);


    if any(params.dp(NNcells) > 0.2)
        params.NN_increased = true;
    else
        params.NN_increased = false;
    end

end
end

%% Core ODE function 

function dy = li(t,y,params)
% y = [d; d_A; n; E; A] concatenated blocks (k each)

global PERIM_STORE T_STORE P_COUNTER

k = params.k;
C_base = params.C_base;
% unpack parameters (local copy)
beta_d = params.beta.d;
beta_n = params.beta.n;
beta_E = params.beta.E;
beta_A = params.beta.A;
alpha_plus = params.alpha.plus;
alpha_minus = params.alpha.minus;
T_s_base = params.T.s;
T_E = params.T.E;
c_s = params.c.s;
c_E = params.c.E;
Kappa_t_base = params.Kappa.t;

% unpack state
y(y<0) = 0;
d   = y(1:k);
d_A = y(k+1:2*k);
n   = y(2*k+1:3*k);
E   = y(3*k+1:4*k);
A   = y(4*k+1:5*k);

%% Determine current perimeter p_t
p0 = params.perimeter0;
p_t = p0;  % default no change

do_delam = isfield(params,'delam') && isstruct(params.delam);
if do_delam
    D = params.delam;
    triggered = false;
    if isfield(D,'trigger') && strcmpi(D.trigger,'time')
        if isfield(D,'time') && t >= D.time
            triggered = true;
        end
    end
    if triggered

        SOP_cells = params.SOP.cells(:);

        if isfield(D,'E_threshold')
            E_active_cells = find(E > D.E_threshold);
        else
            E_active_cells = [];
        end

        % defaults
        if ~isfield(D,'SOP_new_perimeter'), D.SOP_new_perimeter = 0.6; end
        if ~isfield(D,'Ecell_new_perimeter'), D.Ecell_new_perimeter = 1.2; end
        if ~isfield(D,'smooth'), D.smooth = false; end
        if ~isfield(D,'duration'), D.duration = 5; end

        % build target
        p_target = p0;
        if ~isempty(SOP_cells), p_target(SOP_cells) = D.SOP_new_perimeter; end
        if ~isempty(E_active_cells), p_target(E_active_cells) = D.Ecell_new_perimeter.*p0(E_active_cells); end

        if ~D.smooth
            p_t = p_target;
        else
            
            % smooth linear blend
            if strcmpi(D.trigger,'E')
                t0 = t;
            else
                t0 = D.time;
            end

            tau = D.duration;
            if tau <= 0
                p_t = p_target;
            else
                s = min(max((t - t0)/tau, 0), 1);
                p_t = (1-s).*p0 + s.*p_target;
            end
            
        end
    end
end

%% Record perimeter p_t into global recorder
% initialize recorder on first call (or when cleared)
if isempty(P_COUNTER) || P_COUNTER == 0
    P_COUNTER = 1;
end
% initialize arrays if empty
if isempty(PERIM_STORE)
    % start with capacity for 2000 frames (expand if needed)
    initCols = max(2000, 100);
    PERIM_STORE = zeros(k, initCols);
    T_STORE = zeros(1, initCols);
end

% expand if needed
if P_COUNTER > size(PERIM_STORE,2)
    % double capacity
    newCols = size(PERIM_STORE,2) * 2;
    PERIM_STORE(:,end+1:newCols) = 0;
    T_STORE(end+1:newCols) = 0;
end

% store
PERIM_STORE(:,P_COUNTER) = p_t(:);
T_STORE(P_COUNTER) = t;
P_COUNTER = P_COUNTER + 1;

% Recompute contact lengths and normalization using p_t
p = p_t(:);
L = 0.5 * (p + p') .* (C_base > 0);
meanL = mean(L(L>0));

if meanL <= 0
   meanL = 1; 
end
L_norm = L ./ meanL;
C_contact = C_base .* L_norm;

% Recompute cis-inhibition per-cell from p_t (adjusting params.Kappa.InhCombined)
if isfield(params,'cis_size_exponent')
    gamma_cis = params.cis_size_exponent;
else
    gamma_cis = 1;
end

mean_p = mean(p);
cis_factor = (mean_p ./ (p + eps)).^gamma_cis;
cis_factor = cis_factor ./ mean(cis_factor);

if isfield(params,'Kappa') && isfield(params.Kappa,'InhCombined')
    base_kappa = params.Kappa.InhCombined;
else
    base_kappa = 1;
end

if isscalar(base_kappa)
    base_vec = base_kappa .* ones(k,1);
else
    base_vec = base_kappa(:);
    if numel(base_vec)~=k
        base_vec = mean(base_vec).*ones(k,1);
    end
end

Kappa_c = base_vec .* cis_factor;

% Sender-side trans factor (tension) - affects outgoing Delta only
gamma_trans = 1;
if isfield(params,'trans_size_exponent')
    gamma_trans = params.trans_size_exponent;
end
trans_factor = (mean_p ./ (p + eps)).^gamma_trans;
trans_factor = trans_factor ./ mean(trans_factor);
C_sender_weighted = C_contact .* (trans_factor');

% Neighbor contributions
% nneighbor = C_contact * n;              % Notch neighbors (contact-weighted)
nneighbor = C_sender_weighted * n;              % Notch neighbors (contact-weighted)

d_Aneighbor = C_sender_weighted * d_A;  % Activated-Delta from neighbours (sender-weighted)

% Differential equations
dd   = beta_d + alpha_minus .* (d_A) - (alpha_plus .* A + 1 + Kappa_c .* n) .* d;
dd_A = alpha_plus .* A .* d - (alpha_minus + 1 + Kappa_t_base .* nneighbor + Kappa_c .* n) .* d_A;
dn   = beta_n - n - Kappa_t_base .* d_Aneighbor .* n - Kappa_c .* (d + d_A) .* n;

T_s = T_s_base;
signal = (Kappa_t_base ./ T_s) .* d_Aneighbor .* n;

dE = beta_E .* ((signal.^c_s) ./ (1 + signal.^c_s)) - E;

dA = beta_A .* (1 ./ (1 + (E./T_E).^c_E)) - A;

dy = [dd; dd_A; dn; dE; dA];
end
% 
% %% Movie & Perimeter plotting helper
% function F = movielattice_and_perimeter_movies(tout,yout,k,params)
% P = params.P; Q = params.Q;
% 
% % map to real seconds
% real_time_per_model_unit = 150;
% t_real = tout * real_time_per_model_unit;
% 
% % sample every 30 s
% frame_interval = 150;
% t_uniform = 0:frame_interval:max(t_real);          % seconds
% t_frame_model = t_uniform / real_time_per_model_unit;  % model-time unit grid
% 
% % interpolate yout to uniform real-time frames
% y_uniform = interp1(t_real, yout, t_uniform, 'pchip');
% 
% % interpolate perimeter store to same frame times
% if isfield(params,'perimeter_store') && isfield(params,'perimeter_time')
%     perim_time = params.perimeter_time(:);     % column vector
%     perim_store = params.perimeter_store;      % k x N
% 
%     % remove NaN rows 
%     valid = ~isnan(perim_time);
%     perim_time = perim_time(valid);
%     perim_store = perim_store(:,valid);
% 
%     % Find unique times (keep the last occurrence)
%     [perim_time_unique, ia, ~] = unique(perim_time, 'last');
% 
%     % Reorder perimeter store accordingly
%     perim_store_unique = perim_store(:, ia);
% 
%     % interpolate
%     perim_interp = interp1(perim_time_unique, perim_store_unique', ...
%                            t_frame_model, 'pchip', 'extrap');
% else
%     % fallback: constant perimeters
%     perim_interp = repmat(params.perimeter0(:)', length(t_frame_model), 1);
% end
% 
% % prepare video writers
% vE = VideoWriter('SOP_E_movie.avi');
% vE.FrameRate = 1;
% open(vE);
% 
% vP = VideoWriter('SOP_Perimeter_movie.avi');
% vP.FrameRate = 1;
% open(vP);
% 
% % colour scales
% globalMaxE = max(max(yout(:, 3*k+1 : 4*k)));
% SOPcells = params.SOP.cells(:);
% max_E_SOP   = max(y_uniform(:, 3*k + SOPcells));          % E(spl)
% 
% if globalMaxE <= 0, globalMaxE = 1; end
% 
% numFrames = length(t_frame_model);
% 
% for tind = 1:numFrames
%     % E movie frame
%     figE = figure(400); clf('reset'); set(figE,'Visible','off');
%     for i = 1:P
%         for j = 1:Q
%             ind = pq2ind(i,j,P);
%             E_value = y_uniform(tind, 3*k + ind);
%             mycolor_E = min(E_value / globalMaxE, 1);
%             if E_value < 0.015 %max_E_SOP %
%                 color = [1,1,1];
%             else
%                 color = [1 - mycolor_E, 1 - mycolor_E, 1];
%             end
%             plotHexagon(i, j, color);
%             hold on
%         end
%     end
%     axis image; axis off; box off;
%     title(sprintf('E (t = %.0f s)', t_uniform(tind)));
%     frameE = getframe(figE);
%     writeVideo(vE, frameE);
%     close(figE);
% 
%     % Perimeter movie frame
%     figP = figure(401); clf('reset'); set(figP,'Visible','off');
%     % get perimeters for this frame (row vector)
%     perims = perim_interp(tind, :);
%     % normalize to [0,1] for color mapping (white->blue)
%     perim_norm = (perims - min(perim_interp(:))) ./ (max(perim_interp(:)) - min(perim_interp(:)) + eps);
%     for i = 1:P
%         for j = 1:Q
%             ind = pq2ind(i,j,P);
%             pn = perim_norm(ind);
%             colorP = [1 - pn, 1 - pn, 1]; % white -> blueish
%             plotHexagon(i, j, colorP);
%             hold on
%         end
%     end
%     axis image; axis off; box off;
%     title(sprintf('Perimeter (t = %.0f s)', t_uniform(tind)));
%     frameP = getframe(figP);
%     writeVideo(vP, frameP);
%     close(figP);
% end
% 
% close(vE);
% close(vP);
%% Movie & Perimeter plotting helper
function F = movielattice_and_perimeter_movies(tout,yout,k,params)

P = params.P; Q = params.Q;

% map to real seconds
real_time_per_model_unit = 150;
t_real = tout * real_time_per_model_unit;

% sample every 150 s
frame_interval = 150;
t_uniform = 0:frame_interval:max(t_real);          
t_frame_model = t_uniform / real_time_per_model_unit;

% interpolate yout to uniform real-time frames
y_uniform = interp1(t_real, yout, t_uniform, 'pchip');

% interpolate perimeter store to same frame times
if isfield(params,'perimeter_store') && isfield(params,'perimeter_time')
    perim_time = params.perimeter_time(:);
    perim_store = params.perimeter_store;
    
    valid = ~isnan(perim_time);
    perim_time = perim_time(valid);
    perim_store = perim_store(:,valid);
    
    [perim_time_unique, ia, ~] = unique(perim_time, 'last');
    perim_store_unique = perim_store(:, ia);
    
    perim_interp = interp1(perim_time_unique, perim_store_unique', ...
                           t_frame_model, 'pchip', 'extrap');
else
    perim_interp = repmat(params.perimeter0(:)', length(t_frame_model), 1);
end

%% create folders for PDFs
if ~exist('E_frames_pdf','dir')
    mkdir('E_frames_pdf');
end
if ~exist('Perimeter_frames_pdf','dir')
    mkdir('Perimeter_frames_pdf');
end

% prepare video writers
vE = VideoWriter('SOP_E_movie.avi');
vE.FrameRate = 1;
open(vE);

vP = VideoWriter('SOP_Perimeter_movie.avi');
vP.FrameRate = 1;
open(vP);

% colour scales
globalMaxE = max(max(yout(:, 3*k+1 : 4*k)));
SOPcells = params.SOP.cells(:);
max_E_SOP = max(y_uniform(:, 3*k + SOPcells));

if globalMaxE <= 0, globalMaxE = 1; end

numFrames = length(t_frame_model);

for tind = 1:numFrames
    
    %% E MOVIE FRAME
    figE = figure(400); clf('reset'); set(figE,'Visible','off');
    
    for i = 1:P
        for j = 1:Q
            ind = pq2ind(i,j,P);
            E_value = y_uniform(tind, 3*k + ind);
            mycolor_E = min(E_value / globalMaxE, 1);
            
            if E_value < 0.015
                color = [1,1,1];
            else
                color = [1 - mycolor_E, 1 - mycolor_E, 1];
            end
            
            plotHexagon(i, j, color);
            hold on
        end
    end
    
    axis image; axis off; box off;
    title(sprintf('E (t = %.0f s)', t_uniform(tind)));
    % ADD THIS BLOCK
   
    colormap([linspace(1,0,256)', linspace(1,0,256)', ones(256,1)]);
    caxis([0 globalMaxE]);
    cb = colorbar;
    cb.Label.String = 'E value';
    drawnow;

    %% save as PDF
    pdfNameE = fullfile('E_frames_pdf', ...
        sprintf('E_frame_%04d_t_%05ds.pdf', tind, round(t_uniform(tind))));
    exportgraphics(figE, pdfNameE, 'ContentType','vector');

    frameE = getframe(figE);
    writeVideo(vE, frameE);
    close(figE);

    %% PERIMETER MOVIE FRAME
    figP = figure(401); clf('reset'); set(figP,'Visible','off');
    
    perims = perim_interp(tind, :);
    perim_norm = (perims - min(perim_interp(:))) ./ ...
                 (max(perim_interp(:)) - min(perim_interp(:)) + eps);
    
    for i = 1:P
        for j = 1:Q
            ind = pq2ind(i,j,P);
            pn = perim_norm(ind);
            colorP = [1 - pn, 1 - pn, 1];
            plotHexagon(i, j, colorP);
            hold on
        end
    end
    
    axis image; axis off; box off;
    title(sprintf('Perimeter (t = %.0f s)', t_uniform(tind)));
    
    perim_min = min(perim_interp(:));
    perim_max = max(perim_interp(:));
    
    % blue-to-white colormap matching [1 - pn, 1 - pn, 1]
    colormap([linspace(1,0,256)', linspace(1,0,256)', ones(256,1)]);
    
    caxis([perim_min perim_max]);
    cb = colorbar;
    cb.Label.String = 'Perimeter';
    
    drawnow;   % ensures colorbar renders before export

    %% save as PDF
    pdfNameP = fullfile('Perimeter_frames_pdf', ...
        sprintf('Perimeter_frame_%04d_t_%05ds.pdf', tind, round(t_uniform(tind))));
    exportgraphics(figP, pdfNameP, 'ContentType','vector');

    frameP = getframe(figP);
    writeVideo(vP, frameP);
    close(figP);
end

close(vE);
close(vP);
% Perimeter curves mean ± SD for SOPs and nearest neighbours

% Use the interpolated perimeter series per frame (frames x k)
t_minutes = t_frame_model * real_time_per_model_unit / 60;  % minutes

SOPcells = params.SOP.cells(:);
NNcells  = params.SOP.nearest_neighbors(:);
all_cells = 1:k;
others = setdiff(all_cells, [SOPcells; NNcells]);

perim_SOP = perim_interp(:, SOPcells);     % frames x nSOP
perim_NN  = perim_interp(:, NNcells);      % frames x nNN
perim_OTH = perim_interp(:, others);       % frames x nOther

figure(500); clf; hold on;
shadedMeanPlot(t_minutes, perim_SOP, [0.15 0.45 0.7], 'SOP perim');
shadedMeanPlot(t_minutes, perim_NN,  [0.9 0.55 0.0], 'NN perim');
shadedMeanPlot(t_minutes, perim_OTH, [0.6 0.6 0.6], 'Other perim');
xlabel('Time (min)');
ylabel('Perimeter (a.u.)');
title('Perimeter evolution (mean ± SD)');
legend('Location','best');
xlim([min(t_minutes) max(t_minutes)]);
saveas(figure(500),'Perimeters.pdf')

F = 1;
end

%% small helpers
function shadedMeanPlot(t, data, color, labeltxt)
% t: vector (nframes x 1)
% data: nframes x nCells
% color: RGB triple
mu = mean(data, 2);
sd = std(data, 0, 2);
upper = mu + sd;
lower = mu - sd;

% ensure column vectors
t = t(:);
mu = mu(:);
upper = upper(:);
lower = lower(:);

% build patch
xpatch = [t; flipud(t)];
ypatch = [lower; flipud(upper)];
hfill = fill(xpatch, ypatch, color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on;
hline = plot(t, mu, 'Color', color, 'LineWidth', 1.8);
if nargin>=4 && ~isempty(labeltxt)
    set(hline,'DisplayName',labeltxt);
end
end

function clear_global_perimeter_store()
global PERIM_STORE T_STORE P_COUNTER
PERIM_STORE = [];
T_STORE = [];
P_COUNTER = [];
end

%% remaining connectivity helpers and plotHexagon are the same as originals
function C=getconnectivityM(P,Q)
k=P*Q;
C=zeros(k,k);
w=[1/6 1/6 1/6 1/6 1/6 1/6];
for s=1:k
    kneighbour=findneighbourhex(s,P,Q);
    for r=1:6
        C(s,kneighbour(r))=w(r);
    end
end
end

function out = findneighbourhex(ind,P,Q)
[p,q] = ind2pq(ind,P);
out(1) = pq2ind(mod(p,P)+1,q,P);
out(2) = pq2ind(mod(p-2,P)+1,q,P);
qleft = mod(q-2,Q)+1;
qright = mod(q,Q)+1;
if mod(q,2)==1
    pup   = p;
    pdown = mod(p-2,P)+1;
else
    pup   = mod(p,P)+1;
    pdown = p;
end
out(3) = pq2ind(pup,qleft,P);
out(4) = pq2ind(pdown,qleft,P);
out(5) = pq2ind(pup,qright,P);
out(6) = pq2ind(pdown,qright,P);
end

function ind=pq2ind(p,q, P)
ind = p + (q-1)*P;
end

function [p,q]=ind2pq(ind, P)
q = 1+floor((ind-1)/P);
p = ind - (q-1)*P;
end

function plotHexagon(p0,q0,c)
s32 = sqrt(3)/4;
q = q0*3/4;
p = p0*2*s32;
if q0/2 == round(q0/2)
   p = p+s32;
end
x(1) = q-.5; x(2) = q-.25; x(3) = q+.25; x(4) = q+.5; x(5) = q+.25; x(6) = q-.25;
y(1) = p ; y(2) = p+s32; y(3) = p+s32; y(4) = p; y(5) = p-s32; y(6) = p-s32;
c=min(c,ones(1,3));
patch(x,y,c,'linewidth',2);
end
