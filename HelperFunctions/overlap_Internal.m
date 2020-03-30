% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Last Edit:        30.03.2020 (Killian)
%
% Extract the average phase of the component speficied in 'expression' in the domain spedified in 'active_material' from
% the Outersolution specified in 'OuterSolNums' for the soltion specified in 'SolNums'.
% You also need to specify the dataset from where the solution should be taken.
%   std. Values:
%   'expression': 'ewfd.Ex'
%   'active_material': 'OEO'
%   'phases': [-pi/2, -pi/2] [rad]
%   'OuterSolNums': [1, 2]
%   'nr_solution': [1, 1]
%   'dset': 'dset2'
%   'L': 200e-6

function [g0] = overlap_Internal(model, varargin)
%OVERLAP_INTERNAL Summary of this function goes here
%   Detailed explanation goes here
    eps_0 = 8.854e-12; 
    mu_0 = 4*pi*1e-7; 
    c_const = 1/sqrt(eps_0*mu_0); 
    hbar = 1.045718e-34; 

    expr = 'ewfd.Ex';
    active_material = 'OEO';
    phases = [-pi/2; -pi/2];
    OuterSolNums = [1, 2];
    nr_solution = [1, 1];
    dset = 'dset2';
    L = 200e-6; 
    
    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'expression'
                expr = varargin{ii+1};
            case 'active_material'
                active_domain = varargin{ii+1};
                objects = mphgetselection(model.selection(['geom1_' active_domain '_dom']));
                interactive_domain = objects.entities;
            case 'phases'
                phases = varargin{ii+1};
            case 'OuterSolNums'
                OuterSolNums = varargin{ii+1};
            case 'nr_solution'
                nr_solution = varargin{ii+1};
            case 'L'
                L = varargin{ii+1};
            otherwise
        end
    end
    
    if length(OuterSolNums) <= 1
        error('We need two Simulations to calculate this'); 
    end
    
    if length(OuterSolNums) ~= length(nr_solution)
       error('Length of OuterSolNums must be equal to the solutions selected');  
    end
    
    if length(OuterSolNums) > 2
       error('Only two solutions are accepted');  
    end
    
    % Data 1 is the Pump field, Data 2 is the Signal and the idler.
    % Creating a join of both solutions
    model.result().dataset().create("join1", "Join");
    model.result().dataset("join1").set("data", dset);
    model.result().dataset("join1").set("solutions", "one");
    model.result().dataset("join1").set("outersolnum", OuterSolNums(1));
    model.result().dataset("join1").set("solnum", nr_solution(1));
    
    model.result().dataset("join1").set("data2", dset);
    model.result().dataset("join1").set("solutions2", "one");
    model.result().dataset("join1").set("outersolnum2", OuterSolNums(2));
    model.result().dataset("join1").set("solnum2", nr_solution(2));
    model.result().dataset("join1").set("method", "general");
    command = ['exp(i*(', num2str(phases(1)), '))*data1*conj(exp(i*(', num2str(phases(2)), '))*data2)*conj(exp(i*(', num2str(phases(2)), '))*data2)'];
    model.result().dataset("join1").set("expr", command);
    
    % Calculating the overlap integral. Since everywhere except the active
    % domain, the NL Susceptibility should be zero, integrate only there
    overlap = mphint2(model, expr, 'surface', 'dataset', 'join1', 'selection', interactive_domain);
    
    % Calculating the Normalization constants, i.e. Integral of the field
    % over the entire domain should yield exactly one in the end.
    normalization = zeros(3, 1);
    for idx_normalization = 1:length(normalization)-1
       normalization(idx_normalization) = sqrt(mphint2(model, 'ewfd.normE^2', 'surface',...
           'dataset', dset, 'solnum', nr_solution(idx_normalization), ...
           'outersolnum', OuterSolNums(idx_normalization)));
    end
    
    % Assumption of degenerate solution.
    normalization(3) = normalization(2); 
    
    % Calculating the total Normalization constant.
    Atot = 1;
    for idx_normalization = 1:length(normalization)
       Atot = Atot*normalization(idx_normalization);
    end
    
    % Extract effective refractive indices and wavelengths from the
    % simulation
    neff_p = mphglobal(model, 'ewfd.neff', 'dataset', dset, ...
        'outersolnum', OuterSolNums(1), 'solnum', nr_solution(1)); 
    neff_s = mphglobal(model, 'ewfd.neff', 'dataset', dset, ...
        'outersolnum', OuterSolNums(2), 'solnum', nr_solution(2)); 
    neff_i = neff_s; 
    
    wl_p = mphglobal(model, 'ewfd.lambda0', 'dataset', dset, ...
        'outersolnum', OuterSolNums(1), 'solnum', nr_solution(1)); 
    wl_s = mphglobal(model, 'ewfd.lambda0', 'dataset', dset, ...
        'outersolnum', OuterSolNums(2), 'solnum', nr_solution(2)); 
    wl_i = wl_s;
    
    omega_p = 2*pi*c_const/wl_p; 
    omega_s = 2*pi*c_const/wl_s; 
    omega_i = omega_s; 
    B = -1/2*sqrt(omega_p*omega_s*omega_i)/(neff_s*neff_p*neff_i)*sqrt(hbar/eps_0);
    % Nonlinear susceptility coefficient.
    
    chiNL = 1/2*neff_s*neff_i*neff_p^2*0.1e-9;          % 100 pm/V-> r33. Convert to chi=1/2*no^4*r_33
    g0 = B*chiNL/sqrt(L)*overlap/Atot; 
end

