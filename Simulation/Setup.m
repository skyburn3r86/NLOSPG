% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [model] = Setup(varargin)
    %GEOMETRY Setup of Model
    %   Detailed explanation goes here

    % ***** Read Arguments
    % Standard Args
    options = struct(...
        'wSim', {2.5, '[um]'},...
        'hSubstrate', {2, '[um]'},...
        'hWG', {220, '[nm]'},...
        'wWG', {500, '[nm]'},...
        'hOrganic', {200, '[nm]'},...
        'hBuffer', {200, '[nm]'},...
        'hContact', {50, '[nm]'},...
        'r33', {137, '[pm/V]'}, ...
        'lmin', {1100, '[nm]'}, ...
        'lmax', {2300, '[nm]'}, ...
        'Nl', {25, ' '}, ...
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

    % ****** Start
    import com.comsol.model.*
    import com.comsol.model.util.*

    model = ModelUtil.create('Model');
    model.modelPath(pwd);
    model.label('SPGNLO.mph');

    % Set Model Parameters
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        value = options(1).(name);
        unit = options(2).(name);
        if isa(value, 'double')
            value = num2str(value);
        end
        model.param.set(name, [value, ' ', unit]);
    end
    model.component.create('comp1', true);
end

