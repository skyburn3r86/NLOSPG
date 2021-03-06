% Generated through Matlab
% Author:           Christian Haffner
% E-Mail:           christian.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Last Edited:      Christian Haffner
    
function results = calculateDragEffect(model, varargin)
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
expression_field1 = 'ewfd.Ex';
expression_field2 = 'ewfd.Ex';
expression_field3 = 'es.Ex';

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
        case 'title'
            title_str = varargin{ii+1};
        case 'path'
            save_folder = varargin{ii+1};
%         case 'expr_field1'
%             expr_field1 = varargin{ii+1};
%         case 'expr_field2'
%             expr_field2 = varargin{ii+1};
%         case 'expr_field3'
%             expr_field3 = varargin{ii+1};
        otherwise
    end
end

% Caculates the vaccum coupling rate following QuantumConverter Notebook definition based on paper:
% 2018_Superconducting cavity electro-optics A platform for coherent photon conversion between superconducting and photonic circuit
try
    %% Calculating the Lorentz Force
    % TO DO Reduce Expression to single term
    lorentz_forceZ_exp_str = ['1/2*real(-(d(ewfd.Px,x) + d(ewfd.Py,y) + ewfd.Pz*(-1*j*2*pi/ewfd.lambda0*ewfd.neff))*conj(ewfd.Ez)' ...
        '+i*ewfd.omega*(ewfd.Px*conj(ewfd.By) - ewfd.Py*conj(ewfd.Bx)) )'];
    % neglecting the dx = 300pm as integration will counteract it
    % down() lower value at the interface - up() larger value at interface
    lorentz_forceZ_exp_str_surface_yz_right = '1/2*real(-(down(ewfd.Px)-up(ewfd.Px))*conj(ewfd.Ez))';
    lorentz_forceZ_exp_str_surface_yz_left = '1/2*real(-(-down(ewfd.Px)+up(ewfd.Px))*conj(ewfd.Ez))';
    lorentz_forceZ_exp_str_surface_xz_top = '1/2*real(-(down(ewfd.Py)-up(ewfd.Py))*conj(ewfd.Ez))';
    lorentz_forceZ_exp_str_surface_xz_bottom = '1/2*real(-(-down(ewfd.Py)+up(ewfd.Py))*conj(ewfd.Ez))';
    einsteinlaub_force_exp_str = ['1/2*real(ewfd.Px*d(conj(ewfd.Ez),x) + ewfd.Py*d(conj(ewfd.Ez),y)' ....
        '+ ewfd.Pz*(conj(ewfd.Ez)*(1*j*2*pi/ewfd.lambda0*ewfd.neff))'...
        '+i*ewfd.omega*(ewfd.Px*conj(ewfd.By) - ewfd.Py*conj(ewfd.Bx)) )'];
    
    % Saving Force Distribution
    if 1      
        saveSolutionSnapshot(model, 'expression', lorentz_forceZ_exp_str, 'nr_solution', nr_solution,...
            'title', ['Lorentz' title_str], 'path', save_folder);
        saveSolutionSnapshot(model, 'expression', einsteinlaub_force_exp_str, 'nr_solution', nr_solution,...
            'title', ['Einstein' title_str], 'path', save_folder);
    end
    
    % % %
    % Extracting the Surface force for the integration 
    % selection cross-check via
    
    figure(1)
    mphgeom(model, 'geom1', 'Edgelabels', 'on');        
    % TO DO: SELECTIONS
    % Extracting Lorentz Force Surface, 2* for left and right electrode
    lorentz_force_metal_surface_yz_slot = 2*mphint2(model, lorentz_forceZ_exp_str_surface_yz_right, 'line',  'Solnum', nr_solution, 'Selection', [19 21]);
    lorentz_force_metal_surface_yz_outside = 2*mphint2(model, lorentz_forceZ_exp_str_surface_yz_left, 'line', 'Solnum', nr_solution, 'Selection', [14 16]);
    lorentz_force_metal_surface_xz_bottom = 2*mphint2(model, lorentz_forceZ_exp_str_surface_xz_bottom, 'line', 'Solnum', nr_solution, 'Selection', [15]);
    lorentz_force_metal_surface_xz_top = 2*mphint2(model, lorentz_forceZ_exp_str_surface_xz_top, 'line', 'Solnum', nr_solution, 'Selection', [18]);
    % Extracting Lorentz Force
    Lorentz_Force_total = mphint2(model, lorentz_forceZ_exp_str, 'surface', 'Solnum', nr_solution);
    Lorentz_Force_total = Lorentz_Force_total + lorentz_force_metal_surface_yz_slot +lorentz_force_metal_surface_yz_outside + ...
        lorentz_force_metal_surface_xz_bottom + lorentz_force_metal_surface_xz_top;    
    
    Lorentz_Force_metal = mphint2(model, lorentz_forceZ_exp_str, 'surface', 'Solnum', nr_solution, 'Selection', interactive_domain);
    Lorentz_Force_metal = Lorentz_Force_metal + lorentz_force_metal_surface_yz_slot/2 +lorentz_force_metal_surface_yz_outside/2 + ...
        lorentz_force_metal_surface_xz_bottom/2 + lorentz_force_metal_surface_xz_top/2;
    Momentum_change_Lorentz = 2*3e8/n_geff*(Lorentz_Force_total);
  
    Einstein_Laub_Force_total = mphint2(model, einsteinlaub_force_exp_str, 'surface', 'Solnum', nr_solution);
    Einstein_Laub_Force_metal = mphint2(model, einsteinlaub_force_exp_str, 'surface', 'Solnum', nr_solution, 'Selection', interactive_domain);      
    Momentum_change_Einstein_Laub = 2*3e8/n_geff*Einstein_Laub_Force_total;
    
    % detrmining the mode power to nomralize the lorentz field by the
    % power as divPxE is E^2
    Poynting_Vector_time_aver = mphint2(model, 'ewfd.Poavz', 'surface', 'Solnum', nr_solution, 'Selection', 'all');
    dMomentum_dt_Poyinting = 2*2*pi/para.wave0*imag(n_eff)*Poynting_Vector_time_aver;
    
    % Results
    results(1).str = 'F_{Lorentz,total}';
    results(1).unit = '[]'; % TO DO UNITS
    results(1).value = Lorentz_Force_total;          
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'F_{Einstein-Laub,total}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = Einstein_Laub_Force_total;
        
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'dp/dz_{Lorentz}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = Momentum_change_Lorentz;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'dp/dz_{Einstein-Laub}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = Momentum_change_Einstein_Laub;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'dp/dz_{Poynting}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value =  2*2*pi/para.wave0*imag(n_eff(nr_solution))*Poynting_Vector_time_aver;
    
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'Poyinting';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = Poynting_Vector_time_aver;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'F_{Lorentz,metal}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = Lorentz_Force_metal/Poynting_Vector_time_aver;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'F_{Einstein-Laub,metal}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = Einstein_Laub_Force_metal/Poynting_Vector_time_aver;  
    
catch
    error_flag = 1;
end

if error_flag == 1
    
    % Results
    results(1).str = 'F_{Lorentz,total}';
    results(1).unit = '[]'; % TO DO UNITS
    results(1).value = NaN;          
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'F_{Einstein-Laub,total}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = NaN;
        
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'dp/dz_{Lorentz}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = NaN;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'dp/dz_{Einstein-Laub}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = NaN;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'dp/dz_{Poynting}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value =  NaN;
    
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'Poyinting';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = NaN;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'F_{Lorentz,metal}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = NaN;
    last_colmn = size(results,2);
    results(last_colmn+1).str = 'F_{Einstein-Laub,metal}';
    results(last_colmn+1).unit = '[]'; % TO DO UNITS
    results(last_colmn+1).value = NaN;
      
end

end


%% OLD CODE FOR LINE LORENTZ FORCE
% 
%     resolution = 0.15e-9; %[nm]
%     y_cut = para.phot_wg_h + para.coupler_spacing + para.plasm_wg_h/2; % center of plasmonic mim waveguide
%     data_grid(1,:) = [-para.plasm_metal_w-para.plasm_wg_slot_w/2:...
%         resolution:para.plasm_metal_w+para.plasm_wg_slot_w/2];
%     data_grid(2,:) = y_cut;
%  if 0
% %         figure(length(n_eff));
%         model.result('pg2').feature('lngr1').set('expr', lorentz_forceZ_exp_str);
%         model.result('pg2').setIndex('looplevelinput', 'manual', 0);
%         model.result('pg2').setIndex('looplevelindices', nr_solution,0);
%         % Extracting the surface terms
%         P_x = mphinterp(model, 'ewfd.Px', 'coord',data_grid,  'Solnum', nr_solution);
%         E_z = mphinterp(model, 'ewfd.Ez', 'coord',data_grid,  'Solnum', nr_solution);
%         % finding the surface term of the left electrode
%         pos_left = find(abs(data_grid(1,:) - (-para.plasm_wg_slot_w/2 - 1*resolution)) < resolution);
%         pos_right = find(abs(data_grid(1,:) + para.plasm_wg_slot_w/2 - 1*resolution) < resolution);
%         % Calculating the derivatives at the sruface
%         dP_x_dxII = (P_x(max(pos_right)) - P_x(min(pos_left)))/(max(pos_right)-min(pos_left))/resolution;
%         E_z_average = (E_z(max(pos_right)) + E_z(min(pos_left)))/2
%         % comparing with the derivations of Comsol
% %         hold off
% %         dP_x_dx = mpheval(model, 'd(ewfd.Px,x)', 'Solnum', nr_solution, 'Edim', 1, 'Selection', [10 17 22 27 32]);
% %         plot(dP_x_dx.p(1,:),  (dP_x_dx.d1), 'b*');
% %         hold on
% %         plot(-para.plasm_wg_slot_w/2 , real(dP_x_dxII) , '*g');
%         % Calculating the
%         lorentz_force_surface_elec = 1/2*real(-dP_x_dxII*conj(E_z_average));
%     end
%     
%     if 0
%         lorentz_force_line = mpheval(model, lorentz_forceZ_exp_str, 'Solnum', nr_solution, 'Edim', 1, 'Selection', [10 17 22 27 32]);
%         einsteinlaub_force_line = mpheval(model, einsteinlaub_force_exp_str, 'Solnum', nr_solution, 'Edim', 1, 'Selection', [10 17 22 27 32]);
%         figure(3);
%         hold off
%         plot(lorentz_force_line.p(1,:), lorentz_force_line.d1, 'b*')
%         hold on
%         plot(einsteinlaub_force_line.p(1,:), einsteinlaub_force_line.d1, 'r*')
%         plot(-para.plasm_wg_slot_w/2 , lorentz_force_surface_elec , '*g');   
%         xlim([-3*para.plasm_wg_slot_w 3*para.plasm_wg_slot_w])
%         title(['Lorentz n eff:' num2str(round(n_eff(nr_solution).*1000)./1000)])       
%         saveas(gcf,[file.path_model '\Results\' file.str_model_comsol '\' results.string{9} '\' data_str '.jpeg']);       
%     end

