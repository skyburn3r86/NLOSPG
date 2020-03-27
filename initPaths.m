% Generated through Matlab
% Author:           Christian Haffner 
% E-Mail:           christian.haffner1986@gmail.com
% Organization:     ETHZ ITET IEF/NIST
% Date 2020/02/14


function initPaths(modelpath)
%initPaths - loads the path of the base model, functions, class and material data
    global library_path;
    library_path = pwd;
    %% adds the general function not specific to the model 
    folder_content = dir(library_path);
    for jj = 1:length(folder_content)
        if folder_content(jj).isdir == 1 && ~ strcmp(folder_content(jj).name,'.') && ~ strcmp(folder_content(jj).name,'..')
            addpath([library_path '/' folder_content(jj).name]);     
        end
    end
    cd(modelpath)
    
    %% add here model specific folders
    folder_content = dir(modelpath);
    for jj = 1:length(folder_content)
        if folder_content(jj).isdir == 1 && ~ strcmp(folder_content(jj).name,'.') && ~ strcmp(folder_content(jj).name,'..')
            addpath([modelpath '/' folder_content(jj).name]);     
        end
    end 
end
