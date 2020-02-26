% Generated through Matlab
% Author:           Christian Haffner
% E-Mail:           christian.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function results = calculateVaccumCoupling(model, varargin)
%Input Variables: - model = comsol model
% VARARGIN options (extend on need)
% 1. 'type' - defines quantum converter algorithums or spontanous parameteric downconversion
% 2. 'nr_solution' - flag to save or  of the comsol results to plott and save

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
    
    %% scans through varagin. -1 and +1 of for loop due to option/value pairs
    for ii = 1:length(varargin)-1
        switch varargin{ii}
            case 'type'
                type = varargin{ii+1};
            case 'nr_solution'
                nr_solution = varargin{ii+1};
            case 'active_material'
                active_domain = varargin{ii+1};
                objects = mphgetselection(model.selection(['geom1_' active_domain '_dom']));
                interactive_domain = objects.entities;
            otherwise
        end
    end
    
    %% Caculates the vaccum coupling rate following QuantumConverter Notebook definition based on paper: 
    % 2018_Superconducting cavity electro-optics A platform for coherent photon conversion between superconducting and photonic circuit      
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
        '(omega/2*eps0*ewfd.nyy^4*ewfd.Ey*conj(ewfd.Ey)*es.Ey)*r33',...
        'surface', 'solnum', nr_solution, 'selection', interactive_domain);    
    
    % 3. Calcualte the vaccum coupling which includes the 2-omega factor in our definition 
    results.g_0.value = interaction_energy.value/sqrt(normcoeff_optical_energy.value)^2/sqrt(normcoeff_RF.value);
    results.g_0.unit = [interaction_energy.unit{1} '*sqrt(' normcoeff_RF.unit{1} ')*' normcoeff_optical_energy.unit{1}];

    % extracting group refractive index via Sum(Energy)/Sum(PowerFlow) -
    % All-plasmonic Mach-Zehnder Modulator Haffner et al. Nature Photonics
    % (2015)
    [mode_energy.value, mode_energy.unit] = mphint2(model,...
        '(eps0*(ewfd.nxx^2*abs(ewfd.Ex)^2+ewfd.nyy^2*abs(ewfd.Ey)^2+ewfd.nzz^2*abs(ewfd.Ez)^2))',...
        'surface', 'solnum', nr_solution);     
    [mode_power_flow.value, mode_power_flow.unit] = mphint2(model,...
        '(ewfd.Ex*conj(ewfd.Hy)-conj(ewfd.Ey)*(ewfd.Hx))',...
        'surface', 'solnum', nr_solution); 
    
    % imaginary part (loss term) is neglected as 4 orders of magnitude lower
    results.ng.value = real(3e8*mode_energy.value/mode_power_flow.value);    
    results.neff.value = mphglobal(model, 'ewfd.neff', 'dataset', 'dset1', 'outersolnum', 1, 'solnum', nr_solution);
    
    % calculating quality and damping rate including the which includes the 2-omega factor in our definition
    results.Q.value = results.ng.value/(2*abs(imag(results.neff.value)));
    results.gamma.value = mphglobal(model, 'ewfd.omega', 'dataset', 'dset1', 'outersolnum', 1, 'solnum', nr_solution)/results.Q.value;
    results.gamma.unit = 'Hz';
    results.L_prop.value = mphglobal(model, '1/(2*imag(ewfd.neff)*2*pi/wl)', 'dataset', 'dset1', 'outersolnum', 1, 'solnum', nr_solution); 
    results.L_prop.unit = 'm';
end