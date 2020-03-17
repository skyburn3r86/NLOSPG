% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Writes the Data from sim_results into the file with filename with json extension.
%   JSON files have high compatibility to various programming languages, e.g. Java, Python, Matlab. It allows for
%   Structured encoding of Data. However, can get slow for gigantic datasets, as it is not saved in bitcode, but in Ascii.
%   However, for a dataset of under 1000 elements (such as this simulation), the json loads within miliseconds.

function [] = writeToJson(param_list, sim_results, filename)
    toWrite = cell(length(param_list.values), 1);

    for idx_param_list = 1:size(param_list.values, 1)
        element = struct();

        for ii = 1:length(param_list.str)
            element.(param_list.str{ii}) =  struct('value', param_list.values(idx_param_list, ii), 'unit', param_list.unit{ii});
        end

        temp = sim_results{idx_param_list,1};
        for ii = 1:length(temp)
            element.(temp(ii).str) = struct('value', temp(ii).value, 'unit', temp(ii).unit);
        end
        toWrite{idx_param_list} = element;
    end
    if ~contains(filename, 'json')
        filename = [filename '.json'];
    end
    x = jsonencode(toWrite);
    fileID = fopen(filename, 'w');
    fprintf(fileID, x);
    fclose(fileID);
end