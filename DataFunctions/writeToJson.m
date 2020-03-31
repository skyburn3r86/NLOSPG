% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Last Edit:        31.03.2020 (Killian)
%
% Writes the Data from sim_results into the file with filename with json extension. To Parse the parameters of the simulation, use the keyword 'model' to load the model as an argument.
%   JSON files have high compatibility to various programming languages, e.g. Java, Python, Matlab. It allows for
%   Structured encoding of Data. However, can get slow for gigantic datasets, as it is not saved in bitcode, but in Ascii.
%   However, for a dataset of under 1000 elements (such as this simulation), the json loads within miliseconds.

function [] = writeToJson(param_list, sim_results, filename, varargin)
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end
    ArgFlag = false;

    % Check if model was given as argument
    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'model'
                ArgFlag = true;
                model = varargin{ii + 1};
        end
    end
    data = struct();
    params = struct();
    if ArgFlag
        args = mphgetexpressions(model.param);
        for idx_param = 1:size(args, 1)
            params.(args{idx_param, 1}) = struct('value', abs(args{idx_param, 4}), 'phase', angle(args{idx_param, 4}), 'unit', args{idx_param, 5});
        end
    end

    % Parse Simulation Data
    toWrite = cell(length(param_list.values), 1);

    for idx_param_list = 1:size(param_list.values, 1)
        element = struct();

        for ii = 1:length(param_list.str)
            element.(param_list.str{ii}) =  struct('value', param_list.values(idx_param_list, ii), 'unit', param_list.unit{ii});
        end

        idx = ['s', 'p'];           % Signal, Pump. Written to results as array, parse first the first one, then the second one.
        temp = sim_results{idx_param_list,1};
        for ii = 1:length(temp)
            if length(temp(ii).value) > 1
                for kk = 1:length(temp(ii).value)
                    element.([temp(ii).str '_' idx(kk)]) = struct('value', abs(temp(ii).value(kk)), 'phase', angle(temp(ii).value(kk)), 'unit', temp(ii).unit);
                end
            else
                element.(temp(ii).str) = struct('value', abs(temp(ii).value), 'phase', angle(temp(ii).value), 'unit', temp(ii).unit);
            end
        end
        toWrite{idx_param_list} = element;
    end
    if ~contains(filename, 'json')
        filename = [filename '.json'];
    end
    % Write the data to the struct
    data.parameters = params;
    data.results = toWrite;

    % Write to File
    x = jsonencode(data);
    fileID = fopen(filename, 'w');
    fprintf(fileID, x);
    fclose(fileID);
end