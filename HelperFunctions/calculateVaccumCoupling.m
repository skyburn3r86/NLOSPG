% Generated through Matlab
% Author:           Christian Haffner
% E-Mail:           christian.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Last Edited:      Christian Haffner
    
function results = calculateVaccumCoupling(model, varargin)
%Input Variables: - model = comsol model
% VARARGIN options (extend on need)
% 'type' - defines quantum converter algorithums or spontanous parameteric downconversion
% 'nr_solution' - flag to save or  of the comsol results to plott and save
% 'errorFlag' - flag marker 0 (false) or 1 (true) that is handed over from
% previous functions to mark failed simulations. 

% TO DO - add the node finder to check for order of mode
%       - power confinement criteria

    global library_path

    %% checks that paris of input argements are defined
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end
    
    %% default value definition
    type = 'Rf2IR';
    interactive_domain = 1; 
    error_flag = 0;
    
    %% scans through varagin. -1 and +1 of for loop due to option/value pairs
    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'type'
                type = varargin{ii+1};
            case 'nr_solution'
                nr_solution = varargin{ii+1};
            case 'active_material'
                active_domain = varargin{ii+1};
                objects = mphgetselection(model.selection(['geom1_' active_domain '_dom']));
                interactive_domain = objects.entities;
            case 'OuterSolNums'
                OuterSolNums = varargin{ii+1};
            case 'error_flag'
                error_flag = varargin{ii+1};                
            otherwise
        end
    end
    
    if strcmp(type, 'Rf2IR')
        % Caculates the vaccum coupling rate following QuantumConverter Notebook definition based on paper:
        % 2018_Superconducting cavity electro-optics A platform for coherent photon conversion between superconducting and photonic circuit   
        try
            % 1a. Extracting Field normalization constant by calculating the total
            % electrical energy of the optical mode.
            [normcoeff_optical_energy.value, normcoeff_optical_energy.unit] = mphint2(model,...
                '(eps0*(ewfd.nxx^2*(abs(ewfd.Ex)^2 + ewfd.nyy^2*abs(ewfd.Ey)^2 + ewfd.nzz^2*abs(ewfd.Ez)^2)))',...
                'surface', 'solnum', nr_solution);
            % 1b. Nnormalization of RF field is performed via CAPACITANCE -> Energy
            % = 1/2*V_bias^2*C over energy per photon h_bar * omega_rf
            % Note, for spontanous downconversion all three optical fields need to
            % be normalized at various frequencyies. Resulting in an
            % \sqrt(l_cavity) dependency. This is for the RF case included in
            % the capacitance.
            [dummy_RF] = mpheval(model,...
                '(eps0*es.epsilonryy*wWG*lWG/hOEO)*V_bias^2/(hbar*2*pi*f_rf*2)',...
                'solnum', 1, 'dataset', 'dset1', 'selection', interactive_domain);
            normcoeff_RF.value = dummy_RF.d1(1);
            % the unit is 1/m
            normcoeff_RF.unit = [dummy_RF.unit];
            
            % 2. extracts the nonlinear interaction. Note here es.Ey cannot be used
            % as it is part of a different dataset.
            [interaction_energy.value, interaction_energy.unit] = mphint2(model,...
                '(omega/2*eps0*ewfd.nyy^4*ewfd.Ey*(ewfd.Ey)*es.Ey)*r33',...
                'surface', 'solnum', nr_solution, 'selection', interactive_domain);
            
            % 3. Calcualte the vaccum coupling which includes the 2-omega factor in our definition
            results(1).str = 'g_0';
            results(1).unit = '[2\pi Hz]';%[interaction_energy.unit{1} '/sqrt(' normcoeff_RF.unit{1} ')/' normcoeff_optical_energy.unit{1}];
            g_0 = interaction_energy.value/sqrt(normcoeff_optical_energy.value)^2/sqrt(normcoeff_RF.value);
            results(1).value = g_0 ;
            
            
            % extracting group refractive index via Sum(Energy)/Sum(PowerFlow) -
            % Group velocity in lossy periodic structured media - https://journals.aps.org/pra/pdf/10.1103/PhysRevA.82.053825
            % Poyinting Vector density = N/ms = W/m^2 per Area
            % Mode Energy denisty = J/m^3
            [mode_energy.value, mode_energy.unit] = mphint2(model,...
                '(eps0*(ewfd.nxx^2*abs(ewfd.Ex)^2+ewfd.nyy^2*abs(ewfd.Ey)^2+ewfd.nzz^2*abs(ewfd.Ez)^2))',...
                'surface', 'solnum', nr_solution);
            [mode_power_flow.value, mode_power_flow.unit] = mphint2(model,...
                '(ewfd.Ex*conj(ewfd.Hy)-conj(ewfd.Ey)*(ewfd.Hx))',...
                'surface', 'solnum', nr_solution);
            
            % imaginary part (loss term) is neglected as 4 orders of magnitude lower
            ng = real(3e8*mode_energy.value/mode_power_flow.value);
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'ng';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = real(ng);
            
            neff = mphglobal(model, 'ewfd.neff', 'dataset', 'dset1', 'outersolnum', 1, 'solnum', nr_solution);
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'neff';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = neff;
            
            % calculating quality and damping rate including the which includes the 2-omega factor in our definition
            Q = ng/(2*abs(imag(neff)));
            
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'Q';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = Q;
            
            last_colmn = size(results,2);
            gamma = mphglobal(model, 'ewfd.omega', 'dataset', 'dset1', 'outersolnum', 1, 'solnum', nr_solution)/Q;
            results(last_colmn+1).str = 'gamma';
            results(last_colmn+1).unit = '[2\pi x Hz]';
            results(last_colmn+1).value = gamma;
            
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'L_prop_';
            results(last_colmn+1).unit = '[m]';
            results(last_colmn+1).value = mphglobal(model, '1/(2*imag(ewfd.neff)*2*pi/wl)', 'dataset', 'dset1', 'outersolnum', 1, 'solnum', nr_solution);
            
            last_colmn = size(results,2);
            C = 4.*g_0^2/gamma/(2*pi*1e6);
            results(last_colmn+1).str = 'C_reduced';
            results(last_colmn+1).unit = '[for \gamma_{RF} = 2\pi 1MHz]';
            results(last_colmn+1).value = C;
            
        catch
            error_flag = 1;
        end
        
        if error_flag == 1;
            % error during extraction renders simulation invalid --> NaN
            results(1).str = 'g_0';
            results(1).unit = '[2\pi Hz]';
            results(1).value = NaN ;            
            results(last_colmn+1).str = 'ng';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = NaN;            
            results(last_colmn+1).str = 'neff';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = NaN;            
            results(last_colmn+1).str = 'Q';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = NaN;            
            results(last_colmn+1).str = 'gamma';
            results(last_colmn+1).unit = '[2\pi x Hz]';
            results(last_colmn+1).value = NaN;            
            results(last_colmn+1).str = 'L_prop_';
            results(last_colmn+1).unit = '[m]';
            results(last_colmn+1).value = NaN;            
            results(last_colmn+1).str = 'C_reduced';
            results(last_colmn+1).unit = '[for \gamma_{RF} = 2\pi 1MHz]';
            results(last_colmn+1).value = NaN;            
        end

    elseif strcmp(type, 'SPDC')
        % Caculates the vaccum coupling rate following Spontaneous Parametric Downconversion Notebook definition based on paper:
        % Fiorentino, M. et al. (2007). Spontaneous parametric down-conversion in periodically poled KTP waveguides and bulk crystals. Optics Express, 15(12), 7479. 
        % General Quantum Optics theory is described well in Grynberg, G. (2012). Quantization of free radiation. In Introduction to Quantum Optics (pp. 301?324).
        try
            % 1a. Calculating the vacuum coupling rate.
            g0 = overlap_Internal(model, 'nr_solution', nr_solution, 'active_material', active_domain, 'OuterSolNums', OuterSolNums);
            results(1).str = 'g_0';
            results(1).unit = '[2pi Hz]';%[interaction_energy.unit{1} '/sqrt(' normcoeff_RF.unit{1} ')/' normcoeff_optical_energy.unit{1}];
            results(1).value = g0 ;
            
            % Calculating the group refractive inidices, as indicated by
            % chris. 
            % Group velocity in lossy periodic structured media - https://journals.aps.org/pra/pdf/10.1103/PhysRevA.82.053825
            Energy = mphint2(model,...
                '(eps0*(ewfd.nxx^2*abs(ewfd.Ex)^2+ewfd.nyy^2*abs(ewfd.Ey)^2+ewfd.nzz^2*abs(ewfd.Ez)^2))',...
                'surface', 'dataset', 'dset2', 'outersolnum', OuterSolNums(1), 'solnum', nr_solution(1));
            PowerFlow = mphint2(model,...
                '(ewfd.Ex*conj(ewfd.Hy)-conj(ewfd.Ey)*(ewfd.Hx))',...
                'surface', 'dataset', 'dset2', 'outersolnum', OuterSolNums(1), 'solnum', nr_solution(1));
            ng_p = real(3e8*Energy/PowerFlow);
            
            Energy = mphint2(model,...
                '(eps0*(ewfd.nxx^2*abs(ewfd.Ex)^2+ewfd.nyy^2*abs(ewfd.Ey)^2+ewfd.nzz^2*abs(ewfd.Ez)^2))',...
                'surface', 'dataset', 'dset2', 'outersolnum', OuterSolNums(2), 'solnum', nr_solution(2));
            PowerFlow = mphint2(model,...
                '(ewfd.Ex*conj(ewfd.Hy)-conj(ewfd.Ey)*(ewfd.Hx))',...
                'surface', 'dataset', 'dset2', 'outersolnum', OuterSolNums(2), 'solnum', nr_solution(2));
            ng_s = real(3e8*Energy/PowerFlow);
            
            neff_p = mphglobal(model, 'ewfd.neff', 'dataset', 'dset2', ...
                'outersolnum', OuterSolNums(1), 'solnum', nr_solution(1)); 
            neff_s = mphglobal(model, 'ewfd.neff', 'dataset', 'dset2', ...
                'outersolnum', OuterSolNums(2), 'solnum', nr_solution(2)); 
    
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'ng';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = [ng_s, ng_p];
            
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'neff';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = [real(neff_s), real(neff_p)];
            
            % calculating quality and damping rate including the which includes the 2-omega factor in our definition
            Qs = ng_s/(2*abs(imag(neff_s)));
            Qp = ng_p/(2*abs(imag(neff_p)));
            
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'Q';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = [Qs, Qp];
            
            last_colmn = size(results,2);
            gammap = mphglobal(model, 'ewfd.omega', 'dataset', 'dset2', 'outersolnum', 1, 'solnum', nr_solution(1), 'outersolnum', OuterSolNums(1))/Qp;
            gammas = mphglobal(model, 'ewfd.omega', 'dataset', 'dset2', 'outersolnum', 1, 'solnum', nr_solution(2), 'outersolnum', OuterSolNums(2))/Qs;
            results(last_colmn+1).str = 'gamma';
            results(last_colmn+1).unit = '[2\\pi Hz]';
            results(last_colmn+1).value = [gammas, gammap];
            
            last_colmn = size(results,2);
            results(last_colmn+1).str = 'L_prop_';
            results(last_colmn+1).unit = '[m]';
            results(last_colmn+1).value = [-mphglobal(model, '1/(2*imag(ewfd.neff)*2*pi/wl)', 'dataset', 'dset2', 'solnum', nr_solution(2), 'outersolnum', OuterSolNums(2)), ...
                -mphglobal(model, '1/(2*imag(ewfd.neff)*2*pi/wl)', 'dataset', 'dset2', 'solnum', nr_solution(1), 'outersolnum', OuterSolNums(1))];
            
            last_colmn = size(results,2);
            C = 4.*g0^2/gammas;
            results(last_colmn+1).str = 'C_reduced';
            results(last_colmn+1).unit = '[\\gamma_{RF} = 2\\pi 1MHz]';
            results(last_colmn+1).value = C;
            
            last_colmn = size(results,2);
            lc = 1310e-3/(2*abs(real(neff_s) - real(neff_p)));
            results(last_colmn+1).str = 'lc';
            results(last_colmn+1).unit = '[\\mum]';
            results(last_colmn+1).value = lc;
     
        catch
            error_flag = 1;
        end
        
        if error_flag == 1
            % error during extraction renders simulation invalid --> NaN
            results(1).str = 'g_0';
            results(1).unit = '[2\pi Hz]';
            results(1).value = NaN ;            
            results(last_colmn+1).str = 'ng';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = [NaN, NaN];                                
            results(last_colmn+1).str = 'neff';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = [NaN, NaN];            
            results(last_colmn+1).str = 'Q';
            results(last_colmn+1).unit = '[]';
            results(last_colmn+1).value = [NaN, NaN];            
            results(last_colmn+1).str = 'gamma';
            results(last_colmn+1).unit = '[2\pi x Hz]';
            results(last_colmn+1).value = [NaN, NaN];            
            results(last_colmn+1).str = 'L_prop_';
            results(last_colmn+1).unit = '[m]';
            results(last_colmn+1).value = [NaN, NaN];            
            results(last_colmn+1).str = 'C_reduced';
            results(last_colmn+1).unit = '[for \gamma_{RF} = 2\pi 1MHz]';
            results(last_colmn+1).value = NaN; 
            results(last_colmn+1).str = 'lc';
            results(last_colmn+1).unit = '[\\mum]';
            results(last_colmn+1).value = NaN;           
        end
    elseif strcmp(type, 'QED')
        
    else
        error('Type of Coupling rate Calculations not recognized!')        
    end
    
end