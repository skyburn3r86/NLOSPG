% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [nu]=QuantumEvaluation(model, varargin)
    warning('off', 'all')
    options = struct(...
        'wSim', {2.5, '[um]'},...
        'hSubstrate', {2, '[um]'},...
        'hWG', {340, '[nm]'},...
        'wWG', {500, '[nm]'},...
        'hOrganic', {200, '[nm]'},...
        'hBuffer', {200, '[nm]'},...
        'hContact', {50, '[nm]'},...
        'r33', {137, '[pm/V]'}, ...
        'lmin', {1100, '[nm]'}, ...
        'lmax', {2300, '[nm]'}, ...
        'Nl', {25, ' '}, ...
        'wl', {1550, '[nm]'});
    optionNames=fieldnames(options);

    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
       if any(strcmp(inpName,optionNames))
          options(1).(inpName) = pair{2}{1};
          options(2).(inpName) = pair{2}{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end

    for ii = 1:length(optionNames)
        inpName = optionNames{ii};
        if strcmp(options(2).(inpName), '[nm]')
            options(1).(inpName) = options(1).(inpName)*1e-9;
        elseif strcmp(options(2).(inpName), '[um]')
            options(1).(inpName) = options(1).(inpName)*1e-6;
        end
    end

    [x, y, Es, Ei, Ep] = Dressing(model, varargin{:});
    eps0 = 8.854e-12;           % F/m
    mu0 = 4*pi*1e-7;            % H/m
    c_const = 1/sqrt(mu0*eps0); % m/s
    hbar = 1.10545718e-34;      % Joules * s

    chi2_333 = zeros(length(x), length(y));
    chi2 = 1e-9;
    ymin = options.('hSubstrate');
    ymax = ymin + options(1).('hOrganic');

    for ii = 1:length(x)
        for ll = 1:length(y)
            if (x(ii) > 0 && x(ii) < 2.5e-6) && (y(ii) > ymin && y(ii) < ymax)
                chi2_333(ii, ll) = 1;
            end
        end
    end
    Ppump = 1e-3;              % Joules / s
    chi = chi2_333*chi2;
    P = cell(1, 1);
    P2 = cell(1, 1);
    rho = cell(1, 1);
    W = cell(1, 1);
    for ii = 1:length(Es)
        Es1 = Es{ii};
        Ei1 = Ei{ii};
        Ep1 = Ep{ii};
        omegas = Es1.omega;
        neffs = Es1.neff;
        omegai = Ei1.omega;
        neffi = Ei1.neff;
        omegap = Ep1.omega;
        neffp = Ep1.neff;
        l = 10e-3;
        v = 300e6;
        rho{ii} = l^2/(2*pi)^2*real(neffs)*real(neffi)/(hbar*c_const^2);
        [P{ii}, P2{ii}] = SPDC_QS(x, y, chi, Es1, Ei1, Ep1, omegas, neffs, omegai, neffi, omegap, neffp, l, v, Ppump);
        W{ii} = abs(P{ii}).^2*2*pi/hbar*rho{ii}/(hbar*omegas);
        leg{ii} = sprintf('$\\lambda_p$=%.0f, $\\lambda_s$=%.0f, $\\lambda_i$=%.0f', 2*pi*c_const/omegap*1e9, 2*pi*c_const/omegas*1e9, 2*pi*c_const/omegai*1e9);
    end

    N = length(P);
    cmap = autumn(N);
    cmap = cmap(end:-1:1, :);
    maxval = 0;
    nu = zeros(N, 1)
    figure()
    hold on;
    for ii = 1:N
        Wtemp = W{ii};
        Ptemp = P{ii};
        z = linspace(0, l, length(Ptemp));
        plot(z*1e6, abs(Wtemp), 'Color', cmap(ii, :), 'LineWidth', 2);
        if max(abs(Wtemp)) > maxval
           maxval = max(abs(Wtemp));
        end
        nu(ii) = max(abs(Wtemp))/(l*Ppump)
    end
    % make round number out of maxval
    ex = 10^floor(log10(maxval));
    val = ceil(maxval/ex);
    maxval = val*ex;
    hold off;
    grid on;
    grid minor;
    legend(leg, 'Interpreter', 'latex')
    set(gca,'ytick',linspace(0, maxval, 5))
    set(gca,'xtick',linspace(0, max(z)*1e6, 5));
    xlabel('Interaction Length [\mum]', 'FontSize', 14)
    ylabel('Photon Flux Spectral Density [1/s/nm]', 'FontSize', 14);
    saveas(gcf, ['./Figures/' num2str(options(1).('wWG')*1e9) 'nm_' num2str(options(1).('hWG')*1e9) 'nm_QS.png']);
    saveas(gcf, ['./Figures/' num2str(options(1).('wWG')*1e9) 'nm_' num2str(options(1).('hWG')*1e9) 'nm_QS.svg']);
end