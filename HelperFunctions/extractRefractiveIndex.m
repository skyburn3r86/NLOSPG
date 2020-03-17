% Generated through Matlab
% Author:           Christian Haffner
% E-Mail:           christian.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function refractive_index = extractRefractiveIndex(materialName, varargin)
%interpolates the optical propoerties of the file defeind by materialName
% vargin = 0 - extrapolates the wavelength from the current comsol model
% vargin = 1 - takes input value as wavelength value

    if nargin > 2 % 1 and 2 are the mdoel and materialName
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
            error('Arguments needs propertyName/propertyValue pairs')
        end        
        % scans through varagin for wavelength. -1 and +1 of for loop due to
        % option/value pairs
        for ii = 1:length(varargin)-1
            switch varargin{ii}
                case 'wavelength'
                    wavelength = varargin{ii+1};
                case 'model'
                    model = varargin{ii+1};
                    % extract wavelength (wl) from comsol model
                    strWl = char(model.param.get('wl'));
                    wl = str2num(strWl(1:strfind(strWl, '[')-2));
                    prefactor_str = (strWl(strfind(strWl, '['):end));
                    switch prefactor_str
                        case '[pm]'
                            prefactor = 1e-12;
                        case '[nm]'
                            prefactor = 1e-9;
                        case '[um]'
                            prefactor = 1e-6;
                        case '[mm]'
                            prefactor = 1e-3;
                        otherwise
                            prefactor = 1;
                    end
                    wavelength = wl*prefactor;
            end
        end  
    fileID = fopen(materialName,'r');
    dataMaterial = textscan(fileID, '%f%f%f', 'Delimiter', '\t');
    fileID = fclose(fileID);

    if ~isempty(strfind(materialName, 'eps'))
        optical_property = interp1(dataMaterial{1} ,(dataMaterial{2} +1i*dataMaterial{3}),wavelength*1e9, 'linear','extrap');
        refractive_index = sqrt(optical_property);
    elseif ~isempty(strfind(materialName, 'nk'))
        refractive_index = interp1(dataMaterial{1} ,(dataMaterial{2} +1i*dataMaterial{3}),wavelength*1e9, 'linear','extrap');
    else
        refractive_index = -1;
    end

end