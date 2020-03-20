% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Writes the Data from sim_results into the file with filename with json extension.
%   JSON files have high compatibility to various programming languages, e.g. Java, Python, Matlab. It allows for
%   Structured encoding of Data. However, can get slow for gigantic datasets, as it is not saved in bitcode, but in Ascii.
%   However, for a dataset of under 1000 elements (such as this simulation), the json loads within miliseconds.
%   Changelog:
%       v1.0: (18.03.20) Support for Writing Parameters and values into json files
%       v1.1: (20.03.20) Also Writes the non-swept parameters into the file.

function [] = writeToJson(param_list, sim_results, model,  filename)
    Params_NonSwept = mphgetexpressions(model.param);
    Parameters = struct();

    for idx_param_list = 1:size(Params_NonSwept, 1)
        Parameters.(Params_NonSwept{idx_param_list, 1}) = struct('value', abs(Params_NonSwept{idx_param_list, 4}), 'phase', angle(Params_NonSwept{idx_param_list, 4}), 'unit', ['[' Params_NonSwept{idx_param_list, 5} ']']);
    end

    Data = struct();
    toWrite = cell(length(param_list.values), 1);

    for idx_param_list = 1:size(param_list.values, 1)
        element = struct();

        for ii = 1:length(param_list.str)
            element.(param_list.str{ii}) =  struct('value', param_list.values(idx_param_list, ii), 'unit', param_list.unit{ii});
        end
        idx = ['s', 'p']; 
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
    Data.Parameters = Parameters;
    Data.data = toWrite;
    if ~contains(filename, 'json')
        filename = [filename '.json'];
    end
    x = jsonencode(Data);
    fileID = fopen(filename, 'w');
    fprintf(fileID, x);
    fclose(fileID);
end