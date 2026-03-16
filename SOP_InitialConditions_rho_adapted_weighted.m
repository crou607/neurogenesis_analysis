% Script with default initial conditions

% Adapted from Troost et al., 2023


function y0 = SOP_InitialConditions_rho_adapted_weighted(model, params, k)

% Initializes each variable (Delta, d_A, Notch, E, Activator) for every cell.

% INPUTS:
%   model  - string specifying model type (e.g. 'Gaussian-like proneural genes')
%   params - struct of model parameters (includes SOP indices)
%   k      - total number of cells (P*Q)


    % Add small multiplicative noise to all initial conditions 
    noise = 1 + 0.05.*2.*(rand(k,1) - 0.5);  
    % Random factors in [0.95, 1.05] simulate ~5% variation
    % between cells. This breaks perfect symmetry and helps patterning emerge.
    
    % Initialize all state variables to uniform values
    d0   = 1.*ones(k,1);       % Delta
    d_A0 = zeros(k,1);       % Activated Delta
    n0   = 1.*ones(k,1);       % Notch
    E0   = zeros(k,1);         % E(spl) repressor
    A0   = 0.6.*ones(k,1);     % Activator              %%% NOTE %%% 

    % Switch between model types
    switch model
        case 'Mib1 mutual inhibition zone'
            % placeholder (no special initialization in adapted version)
            
        case 'Neur lateral inhibition zone'
            % placeholder (no special initialization in adapted version)
            
        case 'Gaussian-like proneural genes'
            %  Apply Gaussian-like pattern of activator (A0) around SOPs
            Gaussian_steps = [1 0.61 0.14];  
            
            % Assign activator amplitudes
            % params.SOP.cells are defined in SOP_DefaultParams_rho_adapted_weighted.m
            
            % A0(params.SOP.cells)               = Gaussian_steps(1);
            % A0(params.SOP.nearest_neighbors)   = Gaussian_steps(2);
            % A0(params.SOP.next_nearest_neighbors) = Gaussian_steps(3);
           
        otherwise
            warning('Unexpected model type in SOP_InitialConditions_rho_adapted_weighted');
    end

    % Combine into single initial condition vector 
    % Each block of k entries corresponds to one variable:
    % [d, d_A, n, E, A]

    % Multiplying each variable by 'noise' preserves per-cell variation.
    y0 = [d0   .* noise; 
          d_A0 .* noise; 
          n0   .* noise; 
          E0   .* noise; 
          A0   .* noise]; 
    
end


