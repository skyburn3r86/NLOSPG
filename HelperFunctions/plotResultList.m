% Generated by/Last edit by
% Author:           Christian Haffner/Christian Haffner
% E-Mail:           chrisitan.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

%% Funciton: reads out the data(dataYZ_str) from the result list defined by
% paraX.str (and) paraYstr (optional). 
% -> Searches for optimal vlaue in higher dimensions
% If nr_colums > 1(only paraX.str defined)
% If nr_colums > 2(only paraX.str and paraY.str defined)
% list (defined by the rows of parameters_list_values.
% No interpolation possible!!!

function [error_prompt list_vector]= plotResultList(sim_results, para_list, varargin) 

%% input parameters:
% result_List = Matrix (row = nr of simultion, column Y- result y)
% para_list = Matrix (row = nr of  simultionn, column Y- para y)
% VARARGIN options (extend on need)
% 1. 'data_str' - string needs to match results indices
% 2. 'data_FOM' - string of either 'max' or 'min' to search for optimal
% sol. in parameter space
% 3. 'para_str' - cell of strings needs to match results. 
        % 1st element = X-axis,  
        % 2nd element = Y-axis if length(values) > 1  
        % 3rd element = limits by defined values parameter spaced search for optimal solution 
% 4. 'para_values' (Optional) - cell of arrays containg parameters of interest. 
        % Cell length has to match cel length para_str 

% Processing of variable inputs-------------------------------------
if round(length(varargin)/2)~=length(varargin)/2
    error('Arguments needs propertyName/propertyValue pairs')
end
% defeault values
 paper_size = [5 5 6 5.5];
 error_prompt = [];
 legend_flag = 0;


% scans through varagin for wavelength. -1 and +1 of for loop due to
% option/value pairs
for list_loop = 1:2:length(varargin)-1
    if ~iscell(varargin{list_loop})
        switch varargin{list_loop}
            case 'data_str'
                dataYZ.str = varargin{list_loop+1};
                for idx = 1:length(sim_results{1,1})
                    if (strcmp(dataYZ.str, sim_results{1}(idx).str))
                        data_colmn = idx;
                        str_data = sim_results{1}(idx).str;
                        str_unit = sim_results{1}(idx).unit;
                    end
                end
            case 'data_FOM'
                data_FOM.str = varargin{list_loop+1};
            case 'para_str'
                cell_param.str = varargin{list_loop+1};
                for idx_para_str = 1:length( cell_param.str)
                    colm_in_para_list(idx_para_str) = find(strcmp(para_list.str, cell_param.str{idx_para_str}));
                    if isempty(colm_in_para_list(idx_para_str))
                        error(['Call of funciton plotResultList' cell_param.str 'does not match a para_list.str']);
                    else
                        % saves the parameter name
                        para{idx_para_str}.str = cell_param.str{idx_para_str};
                    end
                end
                
            case 'para_unit'
                cell_param.unit = varargin{list_loop+1};
                for idx_para_str = 1:length(cell_param.unit)
                        % saves the parameter name
                        para{idx_para_str}.unit = cell_param.unit{idx_para_str};
                end
            case 'para_value'
                cell_param.values = varargin{list_loop+1};
                for idx_para_str = 1:length( cell_param.str)
                    % uses default parameter value if not defined = empty
                    % cell element
                    if isempty(cell_param.values{idx_para_str})
                        para{idx_para_str}.values =  unique(para_list.values(:,colm_in_para_list(idx_para_str)),'stable');
                    else
                        para{idx_para_str}.values = cell_param.values{idx_para_str};
                    end
                end
            case 'save'
                save_path = varargin{list_loop+1};
            case 'legend_flag'
                legend_flag = varargin{list_loop+1};
            case 'position'
                if isnumeric(varargin{list_loop+1})
                    paper_size = varargin{list_loop+1};
                else
                    error_prompt = ('position needs to be defined by an array of type [posx, posy, width, height] given in centimeters');
                    return
                end
        end
    end
end
clear list_loop
%% Extracting data from list
list_vector = ones(length(para_list.values(:,1)),1);
for idx_parameter = 1:length(para)
    colm_list = colm_in_para_list(idx_parameter);
    row_vector = 0;
    % check if parameter (para) defined by user is within min and maximum of para_list
    if min(para{idx_parameter}.values) > min(para_list.values(:,colm_list))*0.9999
        if max(para{idx_parameter}.values) < max(para_list.values(:,colm_list))*1.0001
            for idx_list = 1:length(para_list.values(:,colm_list))
                rel_error = abs((para{idx_parameter}.values - para_list.values(idx_list,colm_list))./para_list.values(idx_list,colm_list));
                row_vector(idx_list,1) = ~isempty(find(rel_error < 1e-3));
            end
        end
    end
    % contains the row number for the specified parameter set
    list_vector = list_vector & row_vector;
    % error handling
    if sum(row_vector) == 0
        fprintf(['********************************************** \n'...
            'given value:' num2str(para{idx_parameter}.values) ' of paramter:' para{idx_parameter}.str ' does not exist'...
            '********************************************** \n']);
        return
    end
end

clear idx_parameter;
% % consider only values that have been found for all parameters of interest
% list_vector = floor(sum(row_vector)/lenght(para.values)); 
% list_vector  = list_vector(1:4)
counter = 1;
for idx_row_list = 1:length(list_vector)
    if list_vector(idx_row_list)
        result_struct = sim_results{idx_row_list};
        result_struct = result_struct(data_colmn);
        data_plot(counter,1) = result_struct.value;
        para_plot(counter,:) = para_list.values(idx_row_list,:);
        counter = counter + 1;
    end
end
clear counter;
% reshaping data from list into matrix
 for idx_para1 = 1:length(para{1}.values)
        para_1_value = para{1}.values(idx_para1);
        for idx_para2 = 1:length(para{2}.values)
            para_2_value = para{2}.values(idx_para2);
            dummy_value = [];
            for idx_list_loop = 1:length(para_plot)
                current_value_plot_list1 = para_plot(idx_list_loop,colm_in_para_list(1));
                current_value_plot_list2 = para_plot(idx_list_loop,colm_in_para_list(2));
                % search list for parameters 1 and 2 that meet the
                % para_1_values and para_2_values simultanousely
                rel_error_para1 = abs((para_1_value - current_value_plot_list1)./current_value_plot_list1);
                rel_error_para2 = abs((para_2_value -current_value_plot_list2)./current_value_plot_list2);
                if ~isempty(find(rel_error_para1 < 1e-3)) && ~isempty(find(rel_error_para2 < 1e-3))                    
                    if strcmp(data_FOM.str, 'max')
                        if ~isempty(dummy_value)
                            dummy_value = max(data_plot(idx_list_loop), dummy_value);
                        else
                            dummy_value = data_plot(idx_list_loop);
                        end
                    elseif strcmp(data_FOM.str, 'min')
                        if ~isempty(dummy_value)
                            dummy_value = min(data_plot(idx_list_loop), dummy_value);
                        else
                            dummy_value = data_plot(idx_list_loop);
                        end
                    end
                end
            end
        dataZ(idx_para2,idx_para1) = dummy_value;
        end
 end
 
 % 1D plot if second argument has less then 5 values
 if length(para{2}) < 5
     dataX = para{1}.values;
     dataY = dataZ;
     % labeling corresponding to the secon d parameter
        label_str = {};
        for jj = 1:length(para{2}.values)
                label_str{jj} = strrep([ para{2}.str '_' num2str(para{2}.values(jj)) '' para{2}.unit],'_',' ');
        end                
        save_str = [strrep([ str_data ' ' str_unit],'_',' ')];
        save_str = [strrep([ str_data ' ' str_unit],'\',' ')];
        title_str = save_str;
        for jj = 3:length(para)
            if length(para{jj}.values) > 1
                save_str = [save_str '__' para{jj}.str '_[' num2str(min(para{jj}.values)) '-' num2str(max(para{jj}.values)) ']'];
            else
                save_str = [save_str '__' para{jj}.str '_[' num2str(min(para{jj}.values)) ']'];
            end
        end
        
        % generating 1D plot       
        figure('Name', str_data, 'Units', 'centimeters', 'Position', paper_size);
        colors_marker = winter(size(dataY,1)+1);        
        fontsize =9;
        markersize =8;
        linewidth = 1.25;
        colors = lines(2);
        for jj = 1:size(dataY,1)
            plot(dataX, abs(dataY(jj,:))...
                , 'x--', 'MarkerSize',  markersize, 'Linewidth', linewidth, 'Color', colors_marker(jj,:))
            hold on
        end
        xlabel(strrep([ para_list.str{colm_in_para_list(1)} ' ' para_list.unit{colm_in_para_list(1)} ],'_',' '));
        ylabel(strrep([ str_data ' ' str_unit],'_',' '));
        title(title_str);
        if legend_flag
            legend(label_str); 
        end
        grid on;       
        if exist('save_path')
            save_str = [save_path strrep(save_str,'.','pt')];
            saveas(gcf, [save_str '.jpeg']); 
            saveas(gcf, [save_str '.tiff']);
            saveas(gcf, [save_str '.fig']);
        end    
% 2D plot        
else
    dataX = para{1}.values;
    dataY = para{2}.values;
    % generating 1D plot   
    figure('Name', str_data, 'Units', 'centimeters', 'Position', paper_size);
    figure
    surface(dataX, dataY, real(dataZ));
    colorbar;
    save_str = [strrep([ str_data ' ' str_unit],'_',' ')];
    save_str = [strrep([ str_data ' ' str_unit],'\',' ')];
    title_str = save_str;
    for jj = 3:length(para)
        if length(para{jj}.values) > 1
            save_str = [save_str '__' para{jj}.str '_[' num2str(min(para{jj}.values)) '-' num2str(max(para{jj}.values)) ']'];
        else
            save_str = [save_str '__' para{jj}.str '_[' num2str(min(para{jj}.values)) ']'];
        end
    end
    title(title_str);
    xlabel(strrep([ para_list.str{colm_in_para_list(1)} ' ' para_list.unit{colm_in_para_list(1)} ],'_',' '));
    ylabel(strrep([ para_list.str{colm_in_para_list(2)} ' ' para_list.unit{colm_in_para_list(2)} ],'_',' '));
    if exist('save_path')
        save_str = [save_path strrep(save_str,'.','pt')];
        saveas(gcf, [save_str '.jpeg']);  
        saveas(gcf, [save_str '.tiff']); 
        saveas(gcf, [save_str '.fig']); 
    end
end
        
end