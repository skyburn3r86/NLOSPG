% Generated by/Last edit by
% Author:           Christian Haffner/Christian Haffner
% E-Mail:           chrisitan.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

%% Funciton: reads sweep parameter and generates folder for storage of results

function save_path =  generateSavePath(para_sweep)
% input parameters:
% para_sweep = cell of structures containing the parameter sweep information/ or parameter opmtimization limits
% return paramters
% save_path = str continaing the path of the folder
    save_path = [];
    for jj = 1:length(para_sweep)
        if length(para_sweep{jj}.values) > 1
            save_path = [save_path '__' para_sweep{jj}.str];% ...
%                 '_(' num2str(min(para_sweep{jj}.values)) '-' num2str(max(para_sweep{jj}.values)) ')'...
%                 para_sweep{jj}.unit];
        else
            save_path = [save_path '__' para_sweep{jj}.str];% '_(' num2str(min(para_sweep{jj}.values)) ')'...
%                 para_sweep{jj}.unit];
        end
    end
    save_path = strrep(save_path , '.', 'pt');
    
    % generating folder
    [status, msg, msgID] = mkdir(['./Results/' save_path]);
    % check if folder already exist. if true ask user if override should be
    % performed
    if strcmp(msg,'Directory already exists.')
%         selpath = uigetdir('./Results','Folder already exist. Options: Redfine path or overwrite')
        prompt = {'Folder already exist. Options: Redfine path string, overwrite (if name is not changed) or cancel saving:'};
        title = 'File Path';
        dims = [1 length(save_path)+20];
        definput = {save_path};
        save_path = inputdlg(prompt,title,dims,definput); % cancel will return empty str
        [status, msg, msgID] = mkdir(['./Results/' save_path{1}]);
    end
    % transforming cell to string
    save_path = string(save_path);
    save_path = ['./Results/' save_path{1}];
end
