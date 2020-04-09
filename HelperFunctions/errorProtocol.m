% Generated through Matlab
% Author:           Christian Haffner
% E-Mail:           christian.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% last edited:      Christian Haffner
    
function errorProtocol(save_folder, error_at_parameter_str, error_comment_str)
%Input Variables: - model = comsol model
% VARARGIN options (extend on need)
% 1. 'save_folder' - path string where error_protocol.txt can be found. 
% 2. 'error_at_parameter_str' - string of the parameter space for which simulation failed
% 3. 'error_comment_str' - string of detailed comment

 
% open file
    fileID = fopen(['./' save_folder '/error_protocol.txt'], 'a+');
    timestemp = datetime('now');
    fprintf(fileID, [datestr(timestemp) '\n']);
    fprintf(fileID, [error_at_parameter_str '\n']);
    fprintf(fileID, [error_comment_str '\n']);
    fprintf(fileID, ['************************************************\n']);
    fclose(fileID);    
end