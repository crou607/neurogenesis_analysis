% Main script with ODEs and solvers for dynamic perimeter model

% Adapted from Troost et al., 2023

function [yout,tout,params] = SOP_multicell_LI_adapted_DYNAMIC_perimeter_HillsNoRecovery(model,params,y0)

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

% Base connectivity
C_base = getconnectivityM(params.P,params.Q);
params.C_base = C_base;
k = params.P * params.Q;
params.k = k;

% Baseline perimeter
if ~isfield(params,'perimeter') || numel(params.perimeter)~=k
    params.perimeter = ones(k,1);
end
params.perimeter0 = params.perimeter(:);

%% Initial conditions
% --- INITIAL CONDITIONS ---
if nargin < 3 || isempty(y0)
    switch model.ratio
        case 'rho'
            y0_core = SOP_InitialConditions_rho_adapted_weighted(model.zone,params,k);
        case 'epsilon'
            y0_core = SOP_InitialConditions_epsilon_adapted_weighted(model.zone,params,k);
    end
else
    y0_core = y0(:);
end

% Append perimeter if missing
if numel(y0_core) == 5*k
    y0 = [y0_core; params.perimeter0];
elseif numel(y0_core) == 6*k
    y0 = y0_core;
else
    error('y0 has incorrect length: expected %d or %d, got %d',5*k,6*k,numel(y0_core));
end


%% Run solver
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[tout,yout] = ode15s(@li_dynamic,tspan,y0,opts,params);

%% Extract perimeter history
params.perimeter_store = yout(:,5*k+1:6*k)';  % k x Nt
params.perimeter_time  = tout';

%% NN perimeter change detection (unchanged)
params.NN_increased = false;
params.dp = params.perimeter_store(:,end) - params.perimeter_store(:,1);

if isfield(params,'SOP') && isfield(params.SOP,'nearest_neighbors')
    NNcells = params.SOP.nearest_neighbors(:);
    if any(params.dp(NNcells) > 0.2)
        params.NN_increased = true;
    end
end

%% Make movies and plots
F  = movielattice_and_perimeter_movies(tout,yout,k,params);

end

%% CORE ODE 

function dy = li_dynamic(t,y,params)

k = params.k;
C_base = params.C_base;

% unpack state
y(y<0) = 0;
d   = y(1:k);
d_A = y(k+1:2*k);
n   = y(2*k+1:3*k);
E   = y(3*k+1:4*k);
A   = y(4*k+1:5*k);
p   = y(5*k+1:6*k);

p0 = params.perimeter0;

%% PARAMETERS
beta_d = params.beta.d;
beta_n = params.beta.n;
beta_E = params.beta.E;
beta_A = params.beta.A;

alpha_plus = params.alpha.plus;
alpha_minus = params.alpha.minus;

T_s = params.T.s;
T_E = params.T.E;
c_s = params.c.s;
c_E = params.c.E;

Kappa_t = params.Kappa.t;

%% CONTACT GEOMETRY FROM p
L = 0.5 * (p + p') .* (C_base > 0);
meanL = mean(L(L>0));
if meanL <= 0, meanL = 1; end

L_norm = L ./ meanL;
C_contact = C_base .* L_norm;

%% CIS INHIBITION SCALING
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

%% TRANS (SENDER) SCALING
if isfield(params,'trans_size_exponent')
    gamma_trans = params.trans_size_exponent;
else
    gamma_trans = 1;
end

trans_factor = (mean_p ./ (p + eps)).^gamma_trans;
trans_factor = trans_factor ./ mean(trans_factor);

C_sender = C_contact .* (trans_factor');

%% NEIGHBOR INPUTS
% nneighbor   = C_contact * n;

nneighbor   = C_sender * n;

d_Aneighbor = C_sender * d_A;

%% SOP CORE DYNAMICS (UNCHANGED)
dd   = beta_d + alpha_minus .* d_A ...
     - (alpha_plus .* A + 1 + Kappa_c .* n) .* d;

dd_A = alpha_plus .* A .* d ...
     - (alpha_minus + 1 + Kappa_t .* nneighbor + Kappa_c .* n) .* d_A;

dn   = beta_n - n ...
     - Kappa_t .* d_Aneighbor .* n ...
     - Kappa_c .* (d + d_A) .* n;

signal = (Kappa_t ./ T_s) .* d_Aneighbor .* n;

dE = beta_E .* ((signal.^c_s) ./ (1 + signal.^c_s)) - E;
dA = beta_A .* (1 ./ (1 + (E./T_E).^c_E)) - A;


%% PERIMETER DYNAMICS
pmin = 0.05;
pmax = 1.4;

hD = params.perimeter_hill;
hE = params.perimeter_E_hill;

thetaD = params.perimeter_theta_delta;
thetaE = params.perimeter_theta_E;

H_delta = (d_A.^hD) ./ (thetaD^hD + d_A.^hD);
H_E     = (E.^hE)   ./ (thetaE^hE + E.^hE);

dp = ...    
  - params.perimeter_kappa .* H_delta ...
  + H_E;

% hard bounds (numerical safety)
dp(p <= pmin & dp < 0) = 0;
dp(p >= pmax & dp > 0) = 0;
%%
dy = [dd; dd_A; dn; dE; dA; dp];

end

%% HELPER FUNCTIONS (UNCHANGED)

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
x = [q-.5 q-.25 q+.25 q+.5 q+.25 q-.25];
y = [p p+s32 p+s32 p p-s32 p-s32];
c=min(c,ones(1,3));
patch(x,y,c,'linewidth',2);
end

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
%% NEW — create folders for PDFs
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
            
            if E_value < 0.025
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

    % === ADD THIS BLOCK (no change to your colors) ===
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

%% Mean ± SD perimeter plots
t_minutes = t_frame_model * real_time_per_model_unit / 60;

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