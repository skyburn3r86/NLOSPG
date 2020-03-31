% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Last Edit:        31.03.2020 (Killian)
% Reads the content of the json file
%   Detailed Summary Here
function [model_params, sim_results, param_list] = readFromJson(filename)
    if ~contains(filename, '.json')
        error('A JSON file has to be provided')
    end
    fileID = fopen(filename, 'r');
    txt = fscanf(fileID, '%s'); 
    dataStruct = jsondecode(txt);
    fclose(fileID);

    % Done reading
    model_params = dataStruct.parameters;
    data = dataStruct.results;

    sim_results = cell(size(data, 1), 1);
    param_list = struct();
    % Revert to state prior writing to JSON, for postprocessing in MATLAB.
    for idx_data = 1:size(data, 1)              % Parse all results
        % Get result in array element idx_data
        d = data(idx_data)
        % Get the names of this array element
        names = fieldnames(d)
        for idx_names = 1:size(names, 1)                % Parse content of the array element
            tmp_name = names{idx_names};                % Name of the element, e.g. g_0, neff_s
            fieldValues = d.(tmp_name);                 % Value of the element, i.e. value, phase and unit
            tmp_field_names = fieldnames(fieldValues);  % Cell array containing which information about the value is available.

            % The sweep variables are not saved with a phase information, as all are real, however, the data is.
            % Checking for Phase information allows to distinguish between data points and sweep variable informations
            if any(contains(tmp_field_names, 'phase'))
                sim_results{idx_data}.('str') = tmp_name;
                sim_results{idx_data}.('unit') = fieldValues.('unit');
                sim_results{idx_data}.('value') = fieldValues.('value')*exp(1i*fieldValues.('phase'));
            else
                % Check if parameter is already know in the parameter list. If it is not, create an entry
                if ~any(contains('str', fieldnames(param_list)))
                    param_list.str{1, 1} = tmp_name;
                    param_list.unit{1, 1} = fieldValues.('unit');
                end
                if ~any(contains(param_list.str, tmp_name))
                    param_list.str{1, size(param_list.str, 2) + 1} = tmp_name;
                    param_list.unit{1, size(param_list.unit, 2) + 1} = fieldValues.('unit');
                end
                idx_entry = find(strcmp(param_list.str, tmp_name));

                if length(idx_entry) > 1
                    error('IDX_ENTRY is not supposed to be larger than 1. You are not at fault, contact developer of this function.')
                end

                param_list.values(idx_data, idx_entry) = fieldValues.('value');
            end
        end
    end
end