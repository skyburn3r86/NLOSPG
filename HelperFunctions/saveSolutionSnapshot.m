% Generated through Matlab
% Author:           Christian Haffner
% E-Mail:           christian.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function saveSolutionSnapshot(model, varargin)
%Input Variables: - model = comsol model
% VARARGIN options (extend on need)
% 1. 'expression' - expression of the comsol results to plott and save
% 2. 'nr_solution' - flag to save or  of the comsol results to plott and save
% 3. 'title' - string of the figure title will be expanded by neff
% 3. 'path' - string of the path figure will be stored to 

    global library_path

    %% checks that paris of input argements are defined
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end
    
    %% default value definition
    % Threshold value (0.1) needs to be tested for you case.
    expression = 'ewfd.normE';
    nr_solution = 1;
    title_str = 'modeProfile';
    path_str = './Results/FieldProfiles/';
    outerSol = 1; 
    %% scans through varagin. -1 and +1 of for loop due to option/value pairs
        for ii = 1:length(varargin)-1
            switch varargin{ii}
                case 'expression'
                    expression = varargin{ii+1};
                case 'nr_solution'
                    nr_solution = varargin{ii+1};
                case 'title'
                    title_str = varargin{ii+1};
                case 'path'
                    title_str = varargin{ii+1};
                case 'OuterSolNum'
                    outerSol = varargin{ii+1};
                otherwise
            end
        end  
        dset = 'dset1'; 
        if outerSol > 1
           dset = 'dset2';
        end
        %% extract fields
        temp = mpheval(model, expression, 'dataset', dset, 'outersolnum', outerSol, 'solnum', nr_solution);
        neff = mphglobal(model, 'ewfd.neff', 'dataset', dset, 'outersolnum', outerSol, 'solnum', nr_solution);
        title_str = [title_str '_neff_' num2str(real(round(neff*1000)/1000))];
        title_str = strrep(title_str,'.','pt');
        extracted_property = temp.('d1');
        cords = temp.('p');
        x = cords(1, :);
        y = cords(2, :);
        n = 400;
        % define edges
        x_edge=linspace(min(x),max(x),n);
        y_edge=linspace(min(y),max(y),n);
        [X,Y]=meshgrid(x_edge,y_edge);
        indices = [];
        nProp = [];
        
        expression_Z_data = griddata(x,y,extracted_property,X,Y);        
        % Correct NaNs
        II = isnan(expression_Z_data);
        expression_Z_data(II) = 0;
        % Calculate Center of Mass
        dx = X(1, 2) - X(1, 1);
        dy = Y(2, 1) - Y(1, 1);
        
        % Calculate Center of Energy
        mt = sum(sum(expression_Z_data*dx*dy));
        mx = sum(sum(expression_Z_data.*X*dx*dy));
        my = sum(sum(expression_Z_data.*Y*dx*dy));
        xc = mx/mt;
        yc = my/mt;
        
        % Plot and Save the Mode Profile
        figure('Visible', 'off')
        pcolor(X*1e6,Y*1e6,expression_Z_data)
        colormap hot;
        shading interp;
        hold on;
        plot(xc*1e6, yc*1e6, '-k*', 'MarkerSize', 5);
        hold off;
        xlabel('x index [um]');
        ylabel('y index [um]');
        title(strrep(title_str,'_', ' '));
        saveas(gcf, [path_str title_str '.jpeg']);  
        close all
end