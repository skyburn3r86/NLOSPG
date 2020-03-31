% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [model, materials, comsol_parameters] = ModelSetup_Parameters(param_list, para_list_index, varargin)
    %GEOMETRY Setup of Model
    %   Detailed explanation goes here
	global library_path 
    param_list_values = param_list.values;
    param_list_values = param_list_values(para_list_index,:);
    param_list_str = param_list.str;
    
    
	% ***** Read Arguments
    % Standard Args define here the parameters needed to generate model
    % components such as geometry, materials, mesh etc.
    % ASiNx... l0SiNx... Materials defined by analytic expressions
    % e,... c_0 physical constants add new if needed don't delete old
    % wl opertion wavelength of simulation
    % default values are
    comsol_parameters = struct(...
        'wSim', {4, '[um]'},...
        'hSubstrate', {2, '[um]'},...
        'hWG_photonic', {300, '[nm]'},...
        'hWG_metal', {100, '[nm]'},...
        'wWG', {500, '[nm]'},...
        'lWG', {10e-6, '[m]'},...
        'hOrganic', {200, '[nm]'},...
        'hBuffer', {200, '[nm]'},...
        'hOEO', {50, '[nm]'},...
        'hcladding', {2000, '[nm]'},...
        'r33', {137, '[pm/V]'}, ...
        'lmin', {1100, '[nm]'}, ...
        'lmax', {2300, '[nm]'}, ...
        'Nl', {25, ' '}, ...	
        'f_rf', {9e9, '[Hz]'}, ...	
        'n_start', {3, ' '},...
        'nr_modes', {4, ' '},...
        'V_bias', {1, '[V]'}, ...	
        'ASiO2', {0.69617, ' '}, ...
        'BSiO2', {0.40794, ' '}, ...
        'CSiO2', {0.89748, ' '}, ...
        'l0SiO2', {0.068404, '[um]'}, ...
        'l1SiO2', {0.11624, '[um]'}, ...
        'l2SiO2', {9.8961, '[um]'}, ...
        'lp', {0.14586, '[um]'}, ...
        'Gamma', {71428571428571.42, '[Hz]'}, ...
        'ASiNx', {3.0249, ' '}, ...
        'BSiNx', {40314, ' '}, ...
        'l0SiNx', {0.13534, '[um]'}, ...
        'l1SiNx', {1239.842, '[um]'}, ...
        'AAl2O3', {1.2038, '[um]'}, ...
        'BAl2O3', {1.0583, '[um]'}, ...
        'CAl2O3', {5.2807, '[um]'}, ...
        'l0Al2O3', {0.061448, '[um]'}, ...
        'l1Al2O3', {0.1107, '[um]'}, ...
        'l2Al2O3', {17.9266, '[um]'}, ...
        'e', {1.602e-19, '[C]'}, ...
        'hbar', {1.0546e-34, '[J*s]'}, ...
        'eps0', {8.854e-12, '[F/m]'}, ...
        'mu0', {1.2566e-6, '[H/m]'}, ...
        'c_0', {3e8, '[m/s]'}, ...   
        'material_passive_photonic', {3e8, '[m/s]'}, ...    
        'wl', {1550, '[nm]'},...
        'plasmonic_mesh', {150, ' '},...
        'omega', {'2*pi*c_const/wl', ''},...
        'gamma', {'0.001*omega', ''},...
        'sigmaxx', {'-1i*e^2/(4*pi*hbar)*log((-omega + 1i*2*gamma)/(+omega + 1i*2*gamma))', ''},...
        'sigmayy', {0, '[S]'},...
        'sigmazz', {'-1i*e^2/(4*pi*hbar)*log((-omega + 1i*2*gamma)/(+omega + 1i*2*gamma))', ''});	
    optionNames=fieldnames(comsol_parameters);

	% adding new material done by adding line e.g. 'photonic_passive2' , 'filename.xt' 
    % Overlapping geometries -> later entry overrides previous entry (eg. PhotonicWG > Substrate)
    % Note: Plasmonic waveguides materials MUST BE initialized by Metal_xx... 
    %       Electrodes MUST be initialized by Electrodes_xx to enable automating meshing
	materials = struct(...
		'Substrate', 'eps_SiO2_Lemarchand_2013_250nm__2500nm.txt',...
		'PhotonicWG', 'eps_Si3N4_Luke_2015.txt',...
		'OEOWG', 'eps_HD_BB_OH_UniWashington.txt',...
		'OEO', 'eps_HD_BB_OH_UniWashington.txt',...
		'Metal_1', 'eps_Au_IEF_1604.txt',...
		'Metal_2', 'eps_Cu_McPeak.txt',...
		'Graphene', ' ');
    MaterialNames=fieldnames(materials);
	
    % redefine model paramters handed over
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
       if any(strcmp(inpName,optionNames))
          comsol_parameters(1).(inpName) = pair{2}{1};
          comsol_parameters(2).(inpName) = pair{2}{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end  
    
    % stores the current sweep_aparameter as the comsol parameters
    for dimension = 1:length(param_list_str)
       para_name = (param_list_str{dimension});
       if any(strcmp(para_name,optionNames))
          comsol_parameters(1).(para_name) = param_list_values(dimension);
          comsol_parameters(2).(para_name) = param_list.unit{dimension};
       else
          error('%s is not a recognized parameter name',para_name);
       end         
    end

    % ****** Start
    import com.comsol.model.*
    import com.comsol.model.util.*

    model = ModelUtil.create('Model');
    model.modelPath(pwd);
    model.label('OEO_Slot.mph');
    

    % Set Model Parameters
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        str_file = comsol_parameters(1).(name);
        unit = comsol_parameters(2).(name);
        if isa(str_file, 'double')
            str_file = num2str(str_file);
        end
        model.param.set(name, [str_file, ' ', unit]);
    end
    model.component.create('comp1', true);	
    
	% Setup Material - Input wavelength is given in nm input arguments is wl*1e9.
	% I) Analytic (sellmaier, lorentz, drude etc,) defined by model parameters
	model.func.create('Graphene', 'Analytic');
	model.func('Graphene').active(true);
	model.func('Graphene').set('expr', 'imag(sqrt(1-(wl*1e-9)^2/lp^2+1i*Gamma*(wl*1e-9)^3/(2*pi*c_const*lp^2)))');
	model.func('Graphene').set('args', {'wl'});
	model.func('Graphene').set('argunit', 'wl');
	model.func('Graphene').set('plotargs', {'wl' '500' '2000'});
	
	
	% II) Interpolation defined by files in DataIn folder. 3columns with 1st-wavelength, 2nd-real, 3rd-imag 	
	for ii = 1:length(MaterialNames)
        name = MaterialNames{ii};
        str_file = materials.(name); 
        % TO DO NK
        if ~isempty(strfind(materials.(name), '.txt'))
            if ~isempty(strfind(materials.(name), 'eps'))
                % ture if material defined via epsilon
                interpDummy = model.func.create(MaterialNames{ii}, 'Interpolation');
                interpDummy.set('funcs', {['eps' MaterialNames{ii} '_re'] '1'; ['eps' MaterialNames{ii} '_im'] '2'});
                interpDummy.set('source', 'file');
                interpDummy.set('filename', [library_path '\DataIn\' str_file]);
                interpDummy.set('extrap', 'linear');
                interpDummy = [];
            elseif ~isempty(strfind(materials.(name), 'nk'))
                % true if material defined via nk
                interpDummy = model.func.create(MaterialNames{ii}, 'Interpolation');
                interpDummy.set('funcs', {['n' MaterialNames{ii}] '1'; ['k' MaterialNames{ii}] '2'});
                interpDummy.set('source', 'file');
                interpDummy.set('filename', [library_path '\DataIn\' str_file]);
                interpDummy.set('extrap', 'linear');
                interpDummy = [];                
            end
            
        end
    end
end
