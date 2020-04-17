% Generated by/Last edit by
% Author:           Christian Haffner/Christian Haffner
% E-Mail:           chrisitan.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

%% Funciton: reads out ploted field profiles to cross check the simulations
% input: path_name - str of RELATIVE results folder path
% param_list - param_list generated during main_data_acquisition
% list_vector - entries of 1 data used, 0 data not used

function loadFieldProfile(path_name, param_list, list_vector, varargin)

% defeault values
figure_number = 1;

% scans through varagin for wavelength. -1 and +1 of for loop due to
% option/value pairs
for list_loop = 1:2:length(varargin)-1
    if ~iscell(varargin{list_loop})
        switch varargin{list_loop}
            case 'figure_number'
                figure_number = varargin{list_loop+1};
            otherwise
        end
    end
end
%% Extracting data from list
list_entries = find(list_vector ==1);
source_path = cd;
folder = path_name(2:end);
folder = strrep(folder,'/','\');
files = dir([source_path folder]);

% 
f1 = figure(figure_number);
f1.WindowStyle = 'Docked';
f1.Units = 'centimeters';
f1.Position = [5 5 15 15];


for jj = 1:length(list_vector)   
    % Save Comsol Model & Evaluate
    file_str = strrep(param_list.print{list_entries(jj)},'>_<', '');
    file_str = strrep(file_str,'_','');
    file_str = strrep(file_str,'[','');
    file_str = strrep(file_str,']','');
    file_str = strrep(file_str,' ','_');
    for file_idx_folder = 1:length(files)        
        if ~isempty(strfind(files(file_idx_folder).name, file_str))
            image_data = imread([source_path folder files(file_idx_folder).name]);  
            image(image_data);
            disp([file_str '--------- Press any key for next image or ctrl+c to cancel']);
            pause;
        end
    end
end

close(figure_number);
        
end