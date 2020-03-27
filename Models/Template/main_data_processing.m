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
plotResultList(sim_results, param_list,'data_str', 'g_0',...
    'para_str', {'hOEO', 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1550e-9, 300e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);

plotResultList(sim_results, param_list,'data_str', 'Q',...
    'para_str', {para_sweep{1}.str, 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1550e-9, 300e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);

plotResultList(sim_results, param_list,'data_str', 'C_reduced',...
    'para_str', {para_sweep{1}.str, 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1550e-9, 300e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);


%% Plot results
plotResultList(sim_results, param_list,'data_str', 'g_0',...
    'para_str', {'hOEO', 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1650e-9, 300e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);

plotResultList(sim_results, param_list,'data_str', 'gamma',...
    'para_str', {para_sweep{1}.str, 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1650e-9, 300e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);

plotResultList(sim_results, param_list,'data_str', 'C_reduced',...
    'para_str', {para_sweep{1}.str, 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1650e-9, 300e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);

%%
plotResultList(sim_results, param_list,'data_str', 'C_reduced',...
    'para_str', {para_sweep{1}.str, 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1550e-9, 300e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);

% plotResultList(sim_results, param_list,'data_str', 'Q',...
%     'para_str', {para_sweep{1}.str, 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1550e-9, 600e-9},...
%     'data_FOM', 'max', 'save', ['./Results/' save_str '\']);

plotResultList(sim_results, param_list,'data_str', 'C_reduced',...
    'para_str', {para_sweep{1}.str, 'hWG' ,'wl', 'wWG'}, 'para_values', {[], [] ,1550e-9, 900e-9},...
    'data_FOM', 'max', 'save', ['./Results/' save_folder '\']);

%%
plotResultList(sim_results, param_list,'data_str', 'Q', ...
    'para_str', {'wl', para_sweep{1}.str}, 'para_values', {[], []},...
    'data_FOM', 'max', 'save', ['./Results/' workspace_save_str]);