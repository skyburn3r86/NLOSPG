% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [model, refr] = DispersionRelation(model, varargin)
% Function Creates the Parameters for the Dispersion Relation and gives them
% to the Comsol model. Can be given 'plot' as boolean argument to plot the
% dispersion relations
    wl = linspace(1.100, 2.300, 25);
    refr = struct(...
        'Si', zeros(length(wl), 1), ...
        'SiO2', zeros(length(wl), 1), ...
        'Au', zeros(length(wl), 1), ...
        'SiNx', zeros(length(wl), 1), ...
        'Al2O3', zeros(length(wl), 1),...
        'Organics', zeros(length(wl), 1));

    counter = 1;
    options = struct(...
        'plot', false, ...
        'calcRefr', true);
    optionNames=fieldnames(options);

    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
        if any(strcmp(inpName,optionNames))
            options(1).(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
    args = options;
    % Parameters Silicon,
    % Li, H. H. (1980). Refractive index of silicon and germanium and its wavelength and temperature derivatives. Journal of Physical and Chemical Reference Data, 9(3), 561?658. https://doi.org/10.1063/1.555624
    options = struct(...
        'ASi', {3.41696, ' '},...
        'BSi', {0.138497, '[um^2]'},...
        'CSi', {0.013924, '[um^4]'},...
        'DSi', {-0.0000209, '[1/um^2]'},...
        'ESi', {0.000000148, '[1/um^4]'},...
        'l0Si', {0.168, '[um]'});
    optionNames=fieldnames(options);

    for ii = 1:length(optionNames)
        name = optionNames{ii};
        value = options(1).(name);
        unit = options(2).(name);
        if isa(value, 'double')
            value = num2str(value);
        end
        model.param.set(name, [value, ' ', unit]);
    end

    if args.('calcRefr')
        L = 1./(wl.^2 - options(1).('l0Si'));
        y = options(1).('ASi') + options(1).('BSi')*L + options(1).('CSi')*L.^2 + options(1).('DSi').*wl + options(1).('ESi').*wl.^2;
        refr.('Si') = y;
        if args.('plot')
            figure(counter)
            plot(wl, y)
            title('Dispersion Relation Silicon')
            xlabel('Wavelength [um]')
            ylabel('Refractive Index')
            counter = counter + 1;
        end
    end

    % Parameters SiO2
    options = struct(...
        'ASiO2', {0.6961663, ' '},...
        'BSiO2', {0.4079426, ' '},...
        'CSiO2', {0.8974794, ' '},...
        'l0SiO2', {0.0684043, '[um]'},...
        'l1SiO2', {0.1162414, '[um]'},...
        'l2SiO2', {9.8961, '[um]'});
    optionNames=fieldnames(options);
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        value = options(1).(name);
        unit = options(2).(name);
        if isa(value, 'double')
            value = num2str(value);
        end
        model.param.set(name, [value, ' ', unit]);
    end

    if args.('calcRefr')
        y = sqrt(1 + options(1).('ASiO2')*wl.^2./(wl.^2 - options(1).('l0SiO2')^2) + options(1).('BSiO2')*wl.^2./(wl.^2 - options(1).('l1SiO2')^2) + options(1).('CSiO2')*wl.^2./(wl.^2 - options(1).('l2SiO2')^2));
        refr.('SiO2') = y;
        if args.('plot')
            figure(counter)
            plot(wl, y)
            title('Dispersion Relation SiO2')
            xlabel('Wavelength [um]')
            ylabel('Refractive Index')
            counter = counter + 1;
        end
    end

    % Parameters Gold
    % Olmon, R. L., Slovick, B., Johnson, T. W., Shelton, D., Oh, S. H., Boreman, G. D., & Raschke, M. B. (2012). Optical dielectric function of gold. Physical Review B - Condensed Matter and Materials Physics, 86(23), 1?9. https://doi.org/10.1103/PhysRevB.86.235147
    hbar = 6.582119569e-16;
    wp = 8.5/hbar;
    eps0 = 8.85418782e-12;
    mu0 = 4*pi*1e-7;
    c_const = 1/sqrt(eps0*mu0);
    lp = 2*pi*c_const/wp;

    options = struct(...
        'lp', {lp*1e6, ['[um]']},...
        'Gamma', {1/(14e-15), '[Hz]'});

    optionNames=fieldnames(options);
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        value = options(1).(name);
        unit = options(2).(name);
        if isa(value, 'double')
            value = num2str(value);
        end
        model.param.set(name, [value, ' ', unit]);
    end

    if args.('calcRefr')
        y1 = 1 - wl.^2/(options(1).('lp')^2);
        y2 = options(1).('Gamma')*wl.^3/(2*pi*c_const*1e6*options(1).('lp')^2);
        refr.('Au') = sqrt(y1+1i*y2);
        if args.('plot')
            figure(counter)
            plot(wl, real(sqrt(y1+1i*y2)))
            hold on;
            plot(wl, imag(sqrt(y1+1i*y2)))
            hold off;
            title('Dispersion Relation Gold')
            xlabel('Wavelength [um]')
            ylabel('Refractive Index')
            counter = counter + 1;
        end
    end
    % Parameters SiN ! Citation wrong, but in mendeley
    % Luke, K., Yoshitomo, O., Lamont, M. R. E., Gaeta, A. L., & Lipson, M. (2015). Broadband mid-infrared frequency comb generation in a Si3N4 microresonator. Optics Letters, 40(21), 4823?4826. https://doi.org/10.134/OL.40.004823
    options = struct(...
        'ASiNx', {3.0249, ' '},...
        'BSiNx', {40314, ' '},...
        'l0SiNx', {0.1353406, '[um]'},...
        'l1SiNx', {1239.842, '[um]'});
    optionNames=fieldnames(options);
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        value = options(1).(name);
        unit = options(2).(name);
        if isa(value, 'double')
            value = num2str(value);
        end
        model.param.set(name, [value, ' ', unit]);
    end

    if args.('calcRefr')
        y = sqrt(1 + options(1).('ASiNx')*wl.^2./(wl.^2 - options(1).('l0SiNx')^2) + options(1).('BSiNx')*wl.^2./(wl.^2 - options(1).('l1SiNx')^2));
        refr.('SiNx') = y;
        if args.('plot')
            figure(counter)
            plot(wl, y)
            title('Dispersion Relation SiNx')
            xlabel('Wavelength [um]')
            ylabel('Refractive Index')
            counter = counter + 1;
        end
    end

    % Parameters Al2O3 !
    % Malitson, I. H. (1962). Refraction and Dispersion of Synthetic Sapphire. Journal of the Optical Society of America, 52(12), 1377. https://doi.org/10.1364/josa.52.001377
    options = struct(...
        'AAl2O3', {1.203798, ' '},...
        'BAl2O3', {1.058264, ' '},...
        'CAl2O3', {5.280729, ' '},...
        'l0Al2O3', {0.06144821, '[um]'},...
        'l1Al2O3', {0.1106997, '[um]'},...
        'l2Al2O3', {17.92656, '[um]'});
    optionNames=fieldnames(options);
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        value = options(1).(name);
        unit = options(2).(name);
        if isa(value, 'double')
            value = num2str(value);
        end
        model.param.set(name, [value, ' ', unit]);
    end

    if args.('calcRefr')
        y = sqrt(1 + options(1).('AAl2O3')*wl.^2./(wl.^2 - options(1).('l0Al2O3')^2) + options(1).('BAl2O3')*wl.^2./(wl.^2 - options(1).('l1Al2O3')^2) + options(1).('CAl2O3')*wl.^2./(wl.^2 - options(1).('l2Al2O3')^2));
        refr.('Al2O3') = y;
        if args.('plot')
            figure(counter)
            plot(wl, y)
            title('Dispersion Relation Al2O3')
            xlabel('Wavelength [um]')
            ylabel('Refractive Index')
            counter = counter + 1;
        end
    end
    % Parameters Organics - From Chris Material Data
    model.component('comp1').func.create('int1', 'Interpolation');
    model.component('comp1').func('int1').set('source', 'file');
    model.component('comp1').func('int1').setIndex('funcs', 'epsrOrganics', 0, 0);
    model.component('comp1').func('int1').setIndex('funcs', 'epsiOrganics', 1, 0);
    model.component('comp1').func('int1').setIndex('funcs', 1, 0, 1);
    model.component('comp1').func('int1').setIndex('funcs', 2, 1, 1);
    model.component('comp1').func('int1').set('scaledata', 'auto');
    model.component('comp1').func('int1').set('filename', [pwd 'DataIn\eps_HD_BB_OH.txt']);
    model.component('comp1').func('int1').set('extrap', 'linear');
    model.component('comp1').func('int1').set('argunit', 'nm');
    if args.('calcRefr')
        fileID = fopen('DataIn/eps_HD_BB_OH.txt', 'r');
        A = fscanf(fileID, '%f\t%f\t%f\n');
        n = A(2:3:end-1);
        k = A(3:3:end);
        wlI = A(1:3:end-2)*1e-3;
        yy = sqrt(n + 1i*k);
        y = interp1(wlI, yy, wl, 'nearest', 'extrap');
        refr.('Organics') = y;
        if args.('plot')
            figure(counter)
            plot(wlI, real(yy));
            hold on
            plot(wl, real(y));
            plot(wlI, imag(yy));
            plot(wl, imag(y));
            hold off
            title('Dispersion Relation Organics')
            xlabel('Wavelength [um]')
            ylabel('Refractive Index')
            counter = counter + 1;
        end
    end
    % Parameters Graphene
    % Hanson, G. W. (2008). Dyadic Green?s Function for Anisotropic, Non-Local Model of Biased Graphene. IEEE Transactions on Antennas and Propagation, 56(3), 747?757. https://doi.org/10.1109/TAP.2008.917005
    options = struct(...
        'e', {1.602e-19, '[C]'},...
        'hbar', {1.0545718e-34, '[J*s]'},...
        'eps0', {8.854e-12, '[F/m]'},...
        'mu0', {4*pi*1e-7, '[H/m]'},...
        'omega', {'2*pi*c_const/wl', ''},...
        'gamma', {'0.001*omega', ''},...
        'sigmaxx', {'-1i*e^2/(4*pi*hbar)*log((-omega + 1i*2*gamma)/(+omega + 1i*2*gamma))', ''},...
        'sigmayy', {0, '[S]'}, ...
        'sigmazz', {'sigmaxx', ''});
    optionNames=fieldnames(options);
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        value = options(1).(name);
        unit = options(2).(name);
        if isa(value, 'double')
            value = num2str(value);
        end
        model.param.set(name, [value, ' ', unit]);
    end
end