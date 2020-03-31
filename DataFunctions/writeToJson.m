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
    numsols = 2;
    sollabels = cell(numsols, 1);
    for ii = 1:numsols
        sollabels{ii} = num2str(ii); 
    end
    % Check if model was given as argument
    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'model'
                ArgFlag = true;
                model = varargin{ii + 1};
            case 'numsols'
                numsols = varargin{ii + 1};
            case 'sollabels'
                sollabels = varargin{ii + 1}; 
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
        
        temp = sim_results{idx_param_list,1};
        
        % If the str of the value starts with TE or TM, assumes that the
        % struct will be about modes. 
        if strcmp(temp(1).str, 'TE') || strcmp(temp(1).str, 'TM')
            % Is struct, because of differnet Wavelengths
            if numsols > 1
               for idx_sol = 1:numsols
                    for idx_mode = 1:length(temp(1).value{idx_sol})
                        element.([sollabels{idx_sol} temp(1).str '0' num2str(idx_mode - 1)]) = struct('value', abs(temp(1).value{idx_sol}(idx_mode)), 'phase', angle(temp(1).value{idx_sol}(idx_mode)), 'unit', '[a.u]');
                    end
                    
                    if isempty(temp(1).value{idx_sol})
                        element.([sollabels{idx_sol} temp(1).str '00']) = struct('value', 'NaN', 'phase', 'NaN', 'unit', 'NaN'); 
                    end
                    
                    for idx_mode = 1:length(temp(2).value{idx_sol})
                        element.([sollabels{idx_sol} temp(2).str '0' num2str(idx_mode - 1)]) = struct('value', abs(temp(2).value{idx_sol}(idx_mode)), 'phase', angle(temp(2).value{idx_sol}(idx_mode)), 'unit', '[a.u]');
                    end

                    if isempty(temp(2).value)
                        element.([sollabels{idx_sol} temp(2).str '00']) = struct('value', 'NaN', 'phase', 'NaN', 'unit', 'NaN'); 
                    end
               end
            else
                for idx_mode = 1:length(temp(1).value)
                    element.([temp(1).str '0' num2str(idx_mode - 1)]) = struct('value', abs(temp(1).value(idx_mode)), 'phase', angle(temp(1).value(idx_mode)), 'unit', '[a.u]');
                end

                if isempty(temp(1).value)
                    element.([temp(1).str '00']) = struct('value', 'NaN', 'phase', 'NaN', 'unit', 'NaN'); 
                end

                for idx_mode = 1:length(temp(2).value)
                    element.([temp(2).str '0' num2str(idx_mode - 1)]) = struct('value', abs(temp(2).value(idx_mode)), 'phase', angle(temp(2).value(idx_mode)), 'unit', '[a.u]');
                end

                if isempty(temp(2).value)
                    element.([temp(2).str '00']) = struct('value', 'NaN', 'phase', 'NaN', 'unit', 'NaN'); 
                end
            end
        else
            idx = ['s', 'p'];           % Signal, Pump. Written to results as array, parse first the first one, then the second one.

            for ii = 1:length(temp)
                if length(temp(ii).value) == 2
                    for kk = 1:length(temp(ii).value)
                        element.([temp(ii).str '_' idx(kk)]) = struct('value', abs(temp(ii).value(kk)), 'phase', angle(temp(ii).value(kk)), 'unit', temp(ii).unit);
                    end
                else
                    element.(temp(ii).str) = struct('value', abs(temp(ii).value), 'phase', angle(temp(ii).value), 'unit', temp(ii).unit);
                end
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