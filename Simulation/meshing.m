% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = meshing(model, Selection, varargin)
%MESH Meshes the structure. Can be given arguments: 'mat' Si, SiNx and
%'mesh[mat]' for the meshing parameters of the material, 'mesh[mat]x' and
%'mesh[mat]y' for anisotropic meshing (x, y=1 will be isotropic)
    options = struct(...
        'mat', 'Si',...
        'meshAu', {50, '[nm]'},...
        'meshAux', {1, ''},...
        'meshAuy', {5, ''}, ...
        'meshSiNx', {'wl/(2.0*20)', ''},...
        'meshSiNxx', {1, ''},...
        'meshSiNxy', {1, ''}, ...
        'meshSi', {'wl/(3.4*20)', ''},...
        'meshSix', {1, ''},...
        'meshSiy', {1, ''}, ...
        'meshSiO2', {'wl/(1.45*10)', ''},...
        'meshSiO2x', {1, ''},...
        'meshSiO2y', {1, ''}, ...
        'meshAl2O3', {'wl/(1.78*20)', ''},...
        'meshAl2O3x', {1, ''},...
        'meshAl2O3y', {1, ''}, ...
        'meshOrganics', {'wl/(1.78*20)', '[nm]'},...
        'meshOrganicsx', {1, ''},...
        'meshOrganicsy', {1, ''},...,
        'meshAir', {'wl/10', ''},...
        'meshAirx', {1, ''}, ...
        'meshAiry', {1, ''});

    optionNames=fieldnames(options);

    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
       if any(strcmp(inpName,optionNames))
           if isa(pair{2}, 'cell')
                options(1).(inpName) = pair{2}{1};
                options(2).(inpName) = pair{2}{2};
           else
               options(1).(inpName) = pair{2};
               options(2).(inpName) = pair{2};
           end

       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    selectionNames = fieldnames(Selection);
    selectionNames = selectionNames([find(strcmp(selectionNames, 'Au')), find(strcmp(selectionNames, 'Si')),find(strcmp(selectionNames, 'SiNx')), find(strcmp(selectionNames, 'Al2O3')), find(strcmp(selectionNames, 'Organics')), find(strcmp(selectionNames, 'SiO2')), find(strcmp(selectionNames, 'Air'))]);
    counter = 1;
    model.component('comp1').mesh.create('mesh1');
    for ii = 1:length(selectionNames)
        name = selectionNames{ii};
        selection = Selection.(name);
        xscale = options.(['mesh' name 'x']);
        yscale = options.(['mesh' name 'y']);
        if isa(options(1).(['mesh' name]), 'double') || isa(options(1).(['mesh' name]), 'float') || isa(options(1).(['mesh' name]), 'int')
            smax = [num2str(options(1).(['mesh' name])) options(2).(['mesh' name])];
        else
            smax = [options(1).(['mesh' name]) options(2).(['mesh' name])];
        end

        model.component('comp1').mesh('mesh1').create(['ftri', num2str(counter)], 'FreeTri');
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).set('xscale', xscale);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).set('yscale', yscale);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).selection.geom('geom1', 2);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).selection.set(selection);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).create('size1', 'Size');
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).feature('size1').set('hauto', 1);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).feature('size1').set('custom', 'on');
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).feature('size1').set('hmax', smax);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).feature('size1').set('hmaxactive', true);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).feature('size1').set('hmin', 3.86E-11);
        model.component('comp1').mesh('mesh1').feature(['ftri', num2str(counter)]).feature('size1').set('hminactive', false);
        counter = counter + 1;
    end
    model.component('comp1').mesh('mesh1').run;
end

