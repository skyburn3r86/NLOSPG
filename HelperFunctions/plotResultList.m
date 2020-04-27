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

function plotResultList(reduced_results, parameter_axis, reduced_parameters, varargin) 

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
fontsize =9;
markersize =8;
linewidth = 1.25;
colors = lines(2);


% scans through varagin for wavelength. -1 and +1 of for loop due to
% option/value pairs
for list_loop = 1:2:length(varargin)-1
    if ~iscell(varargin{list_loop})
        switch varargin{list_loop}
            case 'data_str'
                data_str_argument = varargin{list_loop+1};
                for idx = 1:length(reduced_results{1,1})
                    if (strcmp(data_str_argument, reduced_results{1}(idx).str))
                        data_idx = idx;
                        data_str = reduced_results{1}(idx).str;
                        unit_str = reduced_results{1}(idx).unit;
                    end
                end
            case 'para_str'
                reduced_parameters_str = varargin{list_loop+1};
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

% generating savepath
save_str = [strrep(data_str,'_',' ')];
% save_str = [strrep([ data_str ' ' unit_str],'\',' ')];
title_str = save_str;
for jj = 1:length(reduced_parameters_str)
    save_str = [save_str '__' reduced_parameters_str{jj}];
end
dt = datestr(now,'yyyy_dd_mm_HHMMSS');
save_str = [save_str '__' dt];

for idx_column = 1:size(reduced_results,2)
    for idx_row = 1:size(reduced_results,1)
        data2plot(idx_row, idx_column) = reduced_results{idx_row, idx_column}(data_idx).value;
    end
end

% 1D plot if second argument has less then 5 values
 if length(parameter_axis{2}.values) < 5
     dataX = parameter_axis{1}.values;
     % labeling corresponding to the secon d parameter
        label_str = {};
        for jj = 1:length(parameter_axis{2}.values)
                label_str{jj} = strrep([ parameter_axis{2}.str '_' num2str(parameter_axis{2}.values(jj)) '' parameter_axis{2}.unit],'_',' ');
        end          
        % generating 1D plot       
        figure('Name', data_str, 'Units', 'centimeters', 'Position', paper_size);
        colors_marker = winter(size(data2plot,1)+1);      
        for jj = 1:size(data2plot,1)
            plot(dataX, (data2plot(jj,:))...
                , 'x--', 'MarkerSize',  markersize, 'Linewidth', linewidth, 'Color', colors_marker(jj,:))
            hold on
        end
        xlabel(strrep([ parameter_axis{1}.str ' ' parameter_axis{1}.unit ],'_',' '));
        ylabel(strrep([ data_str ' ' unit_str],'_',' '));
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
    dataX = parameter_axis{1}.values;
    dataY = parameter_axis{2}.values;
    % generating 1D plot
    figure('Name', data_str, 'Units', 'centimeters', 'Position', paper_size);
    figure
    surface(dataX, dataY, real(data2plot));
    colorbar;
    title_str = save_str;
    title(title_str);
    xlabel(strrep([ parameter_axis{1}.str ' ' parameter_axis{1}.unit ],'_',' '));
    ylabel(strrep([ parameter_axis{2}.str ' ' parameter_axis{2}.unit ],'_',' '));
    if exist('save_path')
        save_str = [save_path strrep(save_str,'.','pt')];
        saveas(gcf, [save_str '.jpeg']);  
        saveas(gcf, [save_str '.tiff']); 
        saveas(gcf, [save_str '.fig']); 
    end
end
        
end