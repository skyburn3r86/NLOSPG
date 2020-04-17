% Generated by/Last edit by
% Author:           Christian Haffner/Christian Haffner
% E-Mail:           chrisitan.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

%% Data Evaluation
close all

% checks if parameters are available in worksepace. if false -> Opens load
% dialog
if ~exist('sim_results') || ~exist('param_list')  || ~exist('sim_results')
    file_str = uigetfile({'*.mat'; '.json'}); 
end

%% highlight available simulation results
% availabe simulation results
available_sim_results = [];
for jj = 1:length(sim_results{1})
    available_sim_results = [available_sim_results ',  ' sim_results{1}(jj).str];
end
fprintf(['************************************************************* \n'...
    'Result Names:' available_sim_results '!\n'...
    '************************************************************* \n']);

% availabe parameters
available_param = [];
for jj = 1:length(param_list.str)
    available_param = [];
    dummy_para_list_values = unique(param_list.values(:,jj));
    for ii = 1:length(dummy_para_list_values)
        available_param = [available_param num2str(dummy_para_list_values(ii))];
        if ii < length(dummy_para_list_values)
            available_param = [available_param ', '];
        end
    end
    available_param = [param_list.str{jj} ': [' available_param ']' param_list.unit{jj}];
fprintf(['Parameter - ' available_param '\n']);
end
fprintf('************************************************************* \n');


%% Plot results
% 1D plot as second para input is < 5 elements
w_slot_array = [];%[2.9333e-08, 4.8667e-08, 6.8e-08, 8.7333e-08, 1.0667e-07, 1.26e-07, 1.4533e-07, 1.6467e-07, 1.84e-07, 2.0333e-07, 2.2267e-07, 2.42e-07, 2.6133e-07, 2.8067e-07, 3e-07];
plotResultList(sim_results, param_list,'data_str', 'g_0',...
    'para_str', {'wSlot', 'hMetal' ,'wl'}, 'para_value', {w_slot_array, [] ,1300e-9},...
    'para_unit', {'[m]', '[m]' , '[m]'},'position', [5 5 16 9],...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\'])

plotResultList(sim_results, param_list,'data_str', 'Q',...
    'para_str', {'wSlot', 'hMetal' ,'wl'}, 'para_value', {w_slot_array, [] ,1300e-9},...
    'para_unit', {'[m]', '[m]' , '[m]'},'position', [5 5 16 9],...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\'])

plotResultList(sim_results, param_list,'data_str', 'C_reduced',...
    'para_str', {'wSlot', 'hMetal' ,'wl'}, 'para_value', {w_slot_array, [] ,1300e-9},...
    'para_unit', {'[m]', '[m]' , '[m]'}, 'position', [5 5 16 9],...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\'])