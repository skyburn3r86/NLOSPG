% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [results] = QuantumEvaluation(model, varargin)
% Function to Evaluate the evolution of the bosonic wavefunction.
%   Takes as arguments the model and the model parameters
    options = struct(...
        'wSim', {2.5, '[um]'},...
        'hSubstrate', {2, '[um]'},...
        'hWG', {340, '[nm]'},...
        'wWG', {500, '[nm]'},...
        'hOrganic', {200, '[nm]'},...
        'hBuffer', {200, '[nm]'},...
        'hContact', {50, '[nm]'},...
        'hAir', {500, '[nm]'},...s
        'r33', {137, '[pm/V]'}, ...
        'lmin', {1100, '[nm]'}, ...
        'lmax', {2300, '[nm]'}, ...
        'Nl', {25, ' '}, ...
        'mat', {0, 'Si'}, ...
        'wl', {[0,0,0], '[nm]'},...
        'zmax', {200, '[um]'}, ...
        'hRidge', {20, '[nm]'}, ...
        'dSlot', {100, '[nm]'}, ...
        'filename', {'', ''}, ...
        'r', {0.5, ''});
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
    % Initialize the propagation distance and vector
    zmin = 0;
    zmax = options.zmax;
    dz = 100e-9;
    zz= linspace(0, zmax, (zmax-zmin)/dz - 1);

    [EpA, EsA, EiA] = ExtractField(model, varargin{:});
    assert((length(EpA) == length(EsA)) && (length(EpA) == length(EiA)), 'Length of the fields are not equal')
    result = struct('lambdas', 0.0, ...
                    'ngs', 0.0, ...
                    'lambdai', 0.0, ...
                    'ngi', 0.0, ...
                    'lambdap', 0.0, ...
                    'ngp', 0.0, ...
                    'g0', 0.0, ...
                    'lps', 0.0, ...
                    'lpi', 0.0, ...
                    'lpp', 0.0, ...
                    'dk', 0.0, ...
                    'l110', 0.0, ...
                    'maxp1', 0.0, ...
                    'l21', 0.0);
    results = cell(length(EpA), 1);
    for nn = 1:length(EpA)
        % Read the first Solution, i.e. Nondegenerate
        Ep = EpA{nn};
        Es = EsA{nn};
        Ei = EiA{nn};

        % Normalize the Fields
        Ep = Ep.normalizeField(zmax);
        Es = Es.normalizeField(zmax);
        Ei = Ei.normalizeField(zmax);

        % Initialize Quantum States
        Es.NPhotons = 10;
        Ei.NPhotons = 10;
        Es = Es.initializeState(1);
        Ei = Ei.initializeState(1);

        % Initialize Pump beam
        tau = dz/Ep.c_const;        % Measure during one propagation slice, also represents the time-bucket size
        Pin = 1e-6;                 % 1 W input Power
        Ep = Ep.setTime(tau);       % Set duration of measurement
        Ep = Ep.setPower(Pin);      % Set Input power

        % Select loss model. Nothing: keeps losses as in model. Lossless: Sets losses to be equal to 0. Lossy: Sets losses to some predefined value.
        lossStr = 'nothing';
        if strcmp(lossStr, 'Lossless')
            chi1is = 0;
            chi1ii = 0;
        elseif strcmp(lossStr, 'Lossy')
            chi1is = -0.0172;
            chi1ii = -0.0220;
        else
            chi1is = Es.chi1i;
            chi1ii = Ei.chi1i;
        end
        Es.chi1i = chi1is;
        Ei.chi1i = chi1ii;

        % Simulate
        [P, Estemp, Eitemp, Eptemp, g0] =  QuantumCalc(Es, Ep, Ei, zz);

        % Plot
        close all;
        record = false;     % Select if result should be written into a video (true) or a 2D Plot (false)
        if record
            writerObj = VideoWriter(['./Videos/' options{1}.('filename') '_M.avi']);
            writerObj.FrameRate = 30;
            open(writerObj);
            myfig = figure(1);
            hold on;

            for ii = 1:length(zz)
                plot(linspace(0, Estemp{ii}.NPhotons, Estemp{ii}.NPhotons + 1), conj(Estemp{ii}.psi).*Estemp{ii}.psi)
                tex = sprintf('z=%.1f', zz(ii)*1e6);
                text(0.8*Es.NPhotons, 0.9, tex)
                p1 = full(conj(Estemp{ii}.psi(2)).*Estemp{ii}.psi(2));
                ylim([0, 1]);
                arrowx = 1;
                arrowy = p1;
                textx = 5;
                texty = 0.5;
                xa = [textx arrowx];
                ya = [texty arrowy];
                [xaf, yaf] = ds2nfu(xa, ya);
                annotation('textarrow', xaf, yaf, 'String', sprintf('P(1) = %.3f', p1));
                xlabel('Photon Number')
                ylabel('Probability')
                drawnow
                myframe = getframe(gca);
                writeVideo(writerObj, myframe);
                clf(myfig);
            end
            hold off;
            close(writerObj);
        end

        if ~record
            [xx,yy] = meshgrid(zz, linspace(1, Estemp{1}.NPhotons-1, Estemp{1}.NPhotons));
            Z = zeros(length(zz), Estemp{1}.NPhotons);
            for ii = 1:length(zz)
                 temp = conj(Estemp{ii}.psi).*Estemp{ii}.psi;
                 Z(ii, :) = temp(2:end);
            end
            figure(1)
            pcolor(xx*1e6, yy, transpose(Z));
            colorbar
            shading interp
            xlabel('Propagation distance z [um]', 'fontsize', 12)
            ylabel('Photon Occupation Number', 'fontsize', 12)

            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',12)
            a = get(gca,'YTickLabel');
            set(gca,'YTickLabel',a,'fontsize',12)
            if ~strcmp(options(1).('filename'), '')
                saveas(gcf, ['./Figures/' options(1).('filename') '_Pr.tiff'], 'tiffn')
            end
        end

        figure(2)
        plot(zz*1e6, P, 'LineWidth', 2, 'Color', [0.8, 0.0, 0.0]);
        xlabel('Propagation distance z [um]', 'fontsize', 12)
        ylabel('Pump Power [W]', 'fontsize', 12)
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12)
        a = get(gca,'YTickLabel');
        set(gca,'YTickLabel',a,'fontsize',12)
        if ~strcmp(options(1).('filename'), '')
            saveas(gcf, ['./Figures/' options(1).('filename') '_P.tiff'], 'tiffn')
        end

        % Get Propagation lengths
        alphas = 2*Es.chi1i/(2*real(Es.neff))*2*pi/Es.wl;
        alphai = 2*Ei.chi1i/(2*real(Ei.neff))*2*pi/Ei.wl;
        alphap = 2*Ep.chi1i/(2*real(Ep.neff))*2*pi/Ep.wl;

        % Get Momentum Mismatch
        dk = Ep.neff*2*pi/Ep.wl - Es.neff*2*pi/Es.wl - Ei.neff*2*pi/Ei.wl;

        % Get Probability Factors
        l110 = zz(min(find(Z(:, 1) > 0.1)));
        maxp1 = max(Z(:, 1));
        l21 = zz(min(find(Z(:, 2) > 0.01))); 

        result.('lps') = 1/(2*alphas);
        result.('lpi') = 1/(2*alphai);
        result.('lpp') = 1/(2*alphap);
        result.('g0') = g0;
        result.('l110') = l110;
        result.('l21') = l21;
        result.('maxp1') = maxp1;
        result.('dk') = dk;

        result.('lambdas') = Es.wl;
        result.('ngs') = Es.wl;
        result.('lambdai') = Ei.wl;
        result.('ngi') = Ei.wl;
        result.('lambdap') = Ep.wl;
        result.('ngp') = Ep.wl;

        results{nn} = result;
end