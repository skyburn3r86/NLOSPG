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
%   'OuterSolNums': 1
%   'nr_solution': 1
%   'dset': 'dset1'
%   'exclude_adj': if part of the active domain should be excluded, specify adjacent domain.
%   (e.g. for slot wg the phase is well defined within the slot, but gets more and more distorted the more one approaches
%   the boundaries, therefore take only the slot as reference for the phase).

function [phase] = PhaseExtraction(model, varargin)
    eps_0 = 8.854e-12;
    mu_0 = 4*pi*1e-7;
    c_const = 1/sqrt(eps_0*mu_0);
    hbar = 1.045718e-34;

    OuterSolNums = 1;
    nr_solution = 1;
    expr = 'ewfd.Ex';
    dset = 'dset1';
    active_domain = 'OEO';
    exclude_adj = '';
    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'nr_solution'
                nr_solution = varargin{ii+1};
            case 'active_material'
                active_domain = varargin{ii+1};
            case 'expression'
                expr = varargin{ii+1};
            case 'OuterSolNums'
                OuterSolNums = varargin{ii+1};
            case 'dset'
                dset = varargin{ii+1};
            case 'exclude_adj'
                exclude_adj = varargin{ii+1}; 
            otherwise
        end
    end
    objects = mphgetselection(model.selection(['geom1_' active_domain '_dom']));
    interactive_domain = objects.entities;

    % Read the Domains to be excluded
    exclusion_domains = []; 
    if ~strcmp(exclude_adj, '')
        if ~isa(exclude_adj, 'cell')
            exclude_adj = {exclude_adj};
        end

        for ii = 1:length(exclude_adj)
            object = mphgetselection(model.selection(['geom1_' exclude_adj{ii} '_dom']));
            domains = object.entities;
            adj_domains = [];
            for idx_dom = 1:length(domains)
                adj_domains = [adj_domains mphgetadj(model, 'geom1', 'domain', 'domain', domains(idx_dom))];
            end
            exclusion_domains = [exclusion_domains, intersect(interactive_domain, adj_domains)];
        end
    end
    interactive_domain = setdiff(interactive_domain, exclusion_domains);

    % Definition of the phase in radiant. arctan(Im(x)/Re(x))
    command = ['atan(imag(', expr, ')/real(', expr, '))'];
    phase = mphmean(model, command, 'surface', 'dataset', dset, 'outersolnum', OuterSolNums, ...
                'solnum', nr_solution, 'selection', interactive_domain);

end