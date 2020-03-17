% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [data] = readFromJson(filename)
% Summary Here
%   Detailed Summary Here
    fileID = fopen(filename, 'r');
    txt = fscanf(fileID, '%s'); 
    data = jsondecode(txt); 
    fclose(fileID);
end