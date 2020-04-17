clear all

% reads the eps files: wavelength, eps_re, eps_im and transform it to nk
% values
formatSpec = '%f %f %f \n';
sizeA = [3 Inf];

% load data
[file_str,path_str] = uigetfile('*.txt');
fileID = fopen([path_str file_str],'r');
x = fscanf(fileID, formatSpec, sizeA);
fclose(fileID);

% transform data
y(1,:) = x(1,:); % wavlength is not changed
y(2,:) = real(sqrt(x(2,:)+1*i*x(3,:))); % n 
y(3,:) = imag(sqrt(x(2,:)+1*i*x(3,:))); % k 

% transfrom file
formatSpec = '%f \t %f \t %f \n';
file_str = strrep(file_str,'eps_','nk_');
fileID = fopen([path_str file_str],'w');
fprintf(fileID, formatSpec, y);
fclose(fileID);