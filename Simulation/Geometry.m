% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model, Selections, BoundarySelection] = Geometry(model, varargin)
%GEOMETRY Generates Geometry of Initialized Parameters (Setup)
    % Get Arguments
    Selections = struct(...
        'SiO2', [], ...
        'SiNx', [], ...
        'Au', [], ...
        'Si', [], ...
        'Organics', [],...
        'Al2O3', []);

    options = struct(...
        'wSim', {2.5, '[um]'},...
        'hSubstrate', {2, '[um]'},...
        'hWG', {340, '[nm]'},...
        'wWG', {500, '[nm]'},...
        'hOrganic', {200, '[nm]'},...
        'hBuffer', {200, '[nm]'},...
        'hContact', {50, '[nm]'},...
        'r33', {137, '[pm/V]'}, ...
        'lmin', {1100, '[nm]'}, ...
        'lmax', {2300, '[nm]'}, ...
        'Nl', {25, ' '}, ...
        'mat', {0, 'Si'}, ...
        'wl', {1550, '[nm]'});
    optionNames=fieldnames(options);

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

    for ii = 1:length(optionNames)
        inpName = optionNames{ii};
        if strcmp(options(2).(inpName), '[nm]')
            options(1).(inpName) = options(1).(inpName)*1e-9;
        elseif strcmp(options(2).(inpName), '[um]')
            options(1).(inpName) = options(1).(inpName)*1e-6;
        end
    end

    % ***************** Create Geometry
    eps = 15e-9;
    model.component('comp1').geom.create('geom1', 2);
    counter = 1;

    % Substrate
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Thermal Oxide');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hSubstrate'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' '0'});
    counter = counter + 1;
    if isempty(Selections.('SiO2'))
        Selections.('SiO2') = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; -eps options(1).('hSubstrate')+eps], 'domain');
    else
        old = Selections.('SiO2');
        new = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; -eps options(1).('hSubstrate')+eps], 'domain');
        Selections.('SiO2') = [old, new];
    end
    % WG
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Waveguide');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wWG' 'hWG'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'(wSim - wWG)/2' 'hSubstrate-hWG'});
    counter = counter + 1;

    % Organics
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Organics');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hOrganic'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate'});
    counter = counter + 1;

    % Buffer
    if ~(options(1).('hBuffer') == 0)
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Buffer Oxide');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hBuffer'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hOrganic'});
        counter = counter + 1;
    end

    % Au
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Gold Contact');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hContact'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hOrganic+hBuffer'});
    counter = counter + 1;

    model.component('comp1').geom('geom1').run;


    if isempty(Selections.(options(2).('mat')))
        Selections.(options(2).('mat')) = mphselectbox(model, 'geom1', [(options(1).('wSim')-options(1).('wWG'))/2-eps (options(1).('wSim')+options(1).('wWG'))/2+eps; options(1).('hSubstrate')-options(1).('hWG')-eps options(1).('hSubstrate')+eps], 'domain');
    else
        old = Selections.(options(2).('mat'));
        new = mphselectbox(model, 'geom1', [(options(1).('wSim')-options(1).('wWG'))/2-eps (options(1).('wSim')+options(1).('wWG'))/2+eps; options(1).('hSubstrate')-options(1).('hWG')-eps options(1).('hSubstrate')+eps], 'domain');
        Selections.(options(2).('mat')) = [old, new];
    end

    if isempty(Selections.('Organics'))
        Selections.('Organics') = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')-eps options(1).('hSubstrate')+options(1).('hOrganic')+eps], 'domain');
    else
        old = Selections.('Organics');
        new = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')-eps options(1).('hSubstrate')+options(1).('hOrganic')+eps], 'domain');
        Selections.('Organics') = [old, new];
    end

    if isempty(Selections.('Al2O3')) && ~(options(1).('hBuffer') == 0)
        Selections.('Al2O3') = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')+options(1).('hOrganic')-eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+eps], 'domain');
    elseif ~(options(1).('hBuffer') == 0)
        old = Selections.('Al2O3');
        new = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')+options(1).('hOrganic')-eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+eps], 'domain');
        Selections.('Al2O3') = [old, new];
    end

    if isempty(Selections.('Au'))
        Selections.('Au') = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')-eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+options(1).('hContact')+eps], 'domain');
    else
        old = Selections.('Au');
        new = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')+options(1).('hOrganic')-eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+eps], 'domain');
        Selections.('Au') = [old, new];
    end

    all = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; -eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+options(1).('hContact')+eps], 'boundary');
    inner = mphselectbox(model, 'geom1', [+eps options(1).('wSim')-eps;eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+options(1).('hContact')-eps], 'boundary');
    BoundarySelection = setdiff(all, inner);
    Graphene = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')-eps options(1).('hSubstrate')+eps], 'boundary');
    TopContact = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')-eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+eps], 'boundary');
    BoundarySelection = struct(...
        'OuterBoundaries', BoundarySelection, ...
        'Graphene', Graphene,...
        'TopContact', TopContact);
end

