% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [model, materials] = ModelSetup_Parameters(varargin)
    %GEOMETRY Setup of Model
    %   Detailed explanation goes here
	global library_path 
    
	% ***** Read Arguments
    % Standard Args define here the parameters needed to generate model
    % components such as geometry, materials, mesh etc.
    % ASiNx... l0SiNx... Materials defined by analytic expressions
    % e,... c_0 physical constants add new if needed don't delete old
    % wl opertion wavelength of simulation
    options = struct(...
        'wSim', {6, '[um]'},...
        'hSubstrate', {4, '[um]'},...
        'hWG_top', {160, '[nm]'},...
        'hWG_bot', {160, '[nm]'},...
        'wWG', {500, '[nm]'},...
        'hOrganic', {200, '[nm]'},...
        'hBuffer', {200, '[nm]'},...
        'hOEO', {50, '[nm]'},...
        'hcladding', {2000, '[nm]'},...
        'r33', {137, '[pm/V]'}, ...
        'lmin', {1100, '[nm]'}, ...
        'lmax', {2300, '[nm]'}, ...
        'Nl', {25, ' '}, ...	
        'f_rf', {9e9, '[Hz]'}, ...	
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
        'plasmonic_mesh', {150, ' '});	
    optionNames=fieldnames(options);
    
    % add golbal parameter structer function to do
    
	
	% adding new material done by adding line e.g. 'photonic_passive2' , 'filename.xt' 
    % Overlapping geometries -> later entry overrides previous entry (eg. PhotonicWG > Substrate)
    % Note: Plasmonic waveguides materials MUST BE initialized by Metal_xx... 
    %       Electrodes MUST be initialized by Electrodes_xx to enable automating meshing
	materials = struct(...
		'Substrate', 'eps_SiO2_Lemarchand_2013_250nm__2500nm.txt',...
		'PhotonicWG', 'eps_Si_Green_2008.txt',...
		'OEO', 'eps_HD_BB_OH_UniWashington.txt',...
		'Metal_1', 'eps_Au_IEF_1604.txt',...
		'Metal_2', 'eps_Cu_McPeak.txt',...
		'Graphene', ' ');
    MaterialNames=fieldnames(materials);
	
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
       if any(strcmp(inpName,optionNames))
          options(1).(inpName) = pair{2}{1};
          options(2).(inpName) = pair{2}{2};
       else
          error('%s is not a recognized parameter name',inpName)
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
        str_file = options(1).(name);
        unit = options(2).(name);
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
        if ~isempty(strfind(materials.(name), '.txt'))
            % ture if interpolation file defined for material
            interpDummy = model.func.create(MaterialNames{ii}, 'Interpolation');
            interpDummy.set('funcs', {['eps' MaterialNames{ii} '_re'] '1'; ['eps' MaterialNames{ii} '_im'] '2'});
            interpDummy.set('funcs', {['eps' MaterialNames{ii} '_re'] '1'; ['eps' MaterialNames{ii} '_im'] '2'});
            interpDummy.set('source', 'file');
            interpDummy.set('filename', [library_path '\' str_file]);
            interpDummy.set('extrap', 'linear');
            interpDummy = [];
        else
            
        end
    end
	
end

