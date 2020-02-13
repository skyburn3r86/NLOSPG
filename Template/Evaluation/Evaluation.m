% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [outerResults] = Evaluation(model, varargin)
%EVALUATION Evaluates the Simulation
%   Detailed explanation goes here
    % Get Data from Arguments
    warning('off', 'all')
    options = struct(...
        'wSim', {2.5, '[um]'},...
        'hSubstrate', {2, '[um]'},...
        'hWG', {340, '[nm]'},...
        'wWG', {500, '[nm]'},...
        'hOrganic', {200, '[nm]'},...
        'hBuffer', {200, '[nm]'},...
        'hContact', {50, '[nm]'},...
        'hAir', {500, '[nm]'},...
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
    outerResults = cell(options(1).('Nl'), 1);
    dl = (options(1).('lmax')-options(1).('lmin'))/(options(1).('Nl')-1);

    % Evaluate Data, Sorting and Mode Selection
    for iter = 1:options(1).('Nl')
        options(1).('wl') = options(1).('lmin') + dl*(iter-1);

        % Get the effective refractive indeces
        temp = mpheval(model, 'ewfd.neff', 'dataset', 'dset3', 'outersolnum', iter);
        neff = temp.('d1');
        neff = neff(:, 1);
        [~, sI] = sort(imag(neff));
        neff = neff(sI);

        temp = mpheval(model, 'ewfd.normE', 'dataset', 'dset3', 'outersolnum', iter);
        normE = temp.('d1');
        cords = temp.('p');
        x = cords(1, :);
        y = cords(2, :);
        normE = normE(sI, :);

        temp = mpheval(model, 'ewfd.nxx', 'dataset', 'dset3', 'outersolnum', iter);
        nRefr = temp.('d1');
        nRefr = nRefr(sI, :);

        n = 200;

        x_edge=linspace(min(x),max(x),n);
        y_edge=linspace(min(y),max(y),n);
        [X,Y]=meshgrid(x_edge,y_edge);
        indices = [];
        nProp = [];
        for ii = 1:length(sI)
            Z = griddata(x,y,normE(ii, :),X,Y);
            Zn = griddata(x, y, nRefr(ii, :), X, Y);

            % Correct NaNs
            II = isnan(Z);
            Z(II) = 0;
            II = isnan(Zn);
            Zn(II) = 0;
            % Calculate Center of Mass
            dx = X(1, 2) - X(1, 1);
            dy = Y(2, 1) - Y(1, 1);

            % Calculate Center of Energy
            mt = sum(sum(Z*dx*dy));
            mx = sum(sum(Z.*X*dx*dy));
            my = sum(sum(Z.*Y*dx*dy));

            xc = mx/mt;
            yc = my/mt;
            if ~((xc > (options(1).('wSim') - options(1).('wWG'))/2) && (xc < (options(1).('wSim') + options(1).('wWG'))/2) && ...
                (yc > (options(1).('hSubstrate') - options(1).('hWG'))) && (yc < (options(1).('hSubstrate') + options(1).('hOrganic'))))
                continue;
            end

            % Select array elements to be evaluated (inside WG), ratio is the Guided Power
            I = (X > (options(1).('wSim') - options(1).('wWG'))/2) & (X < (options(1).('wSim') + options(1).('wWG'))/2) & ...
                (Y > (options(1).('hSubstrate') - options(1).('hWG'))) & (Y < (options(1).('hSubstrate')));
            Ptot = abs(sum(sum(sqrt(Zn).*Z.^2)));
            Pinner = abs(sum(sum(sqrt(Zn(I)).*Z(I).^2)));
            ratio = Pinner/Ptot;

            % Select Modes
            n1 = max(real(Zn(I)));
            n2 = 1.454;
            f = (n1^2 - n2^2)/(2*n1);
            r = f*options(1).('wWG')*options(1).('hWG')*4/(pi*options(1).('wl')^2);
            xi = (x_edge > (options(1).('wSim') - options(1).('wWG'))/2 + 5*dx) & (x_edge < (options(1).('wSim') + options(1).('wWG'))/2 - 5*dx);
            yi = find((y_edge > (options(1).('hSubstrate') - options(1).('hWG')/2)) & (y_edge < (options(1).('hSubstrate') - options(1).('hWG')/2) + dy));
            if length(yi) > 1
                yi = yi(1);
            end

            pk = findpeaks(Z(180, :));      % Local Maximas underneath the Metal
            pk2 = findpeaks(Z(yi, xi));     % Local Maximas in the Waveguide
            if ~((ratio > r && length(pk2) == 1) || (length(pk) == 1 && length(pk2) == 1))
                continue
            end

            % Plot and Save the Mode Profile
            figure('Visible', 'off')
            pcolor(X*1e6,Y*1e6,Z)
            colormap hot
            shading interp
            hold on;
            plot(xc*1e6, yc*1e6, '-k*', 'MarkerSize', 10)
            hold off;
            xlabel('x index')
            ylabel('y index')
            title([num2str(options(1).('wl')) 'nm: E for ' num2str(real(neff(ii))) ' + ' num2str(imag(neff(ii))) +'i'])
            saveas(gcf, ['./Figures/FieldProfiles/' num2str(options(1).('wWG')*1e9) 'nm_' num2str(options(1).('hWG')*1e9) 'nm_' num2str(options(1).('hOrganic')*1e9) 'nm_' num2str(options(1).('wl')*1e9) 'nm_' num2str(real(neff(ii)), 4) '-i' num2str(abs(imag(neff(ii))), 4) '_Enorm.png'])
            indices = [indices sI(ii)];
            nProp = [nProp neff(ii)];
        
        end

        % Get Fields for Polarization Determination
        pol = cell(length(indices), 1);
        temp = mpheval(model, 'ewfd.Ex', 'dataset', 'dset3', 'outersolnum', iter);
        Ex = temp.('d1');
        Ex = Ex(indices, :);

        temp = mpheval(model, 'ewfd.Ey', 'dataset', 'dset3', 'outersolnum', iter);
        Ey = temp.('d1');
        Ey = Ey(indices, :);

        temp = mpheval(model, 'ewfd.normE', 'dataset', 'dset3', 'outersolnum', iter);
        normE = temp.('d1');
        normE = normE(indices, :);

        for ii = 1:length(indices)
            Zx=griddata(x,y, Ex(ii, :), X, Y);
            II = isnan(Zx);
            Zx(II) = 0;

            Zy=griddata(x,y, Ey(ii, :), X, Y);
            II = isnan(Zy);
            Zy(II) = 0;

            Zt=griddata(x,y, normE(ii, :), X, Y);
            II = isnan(Zt);
            Zt(II) = 0;

            % Calculate TE Fraction -- Following Lumerical Definition
            yy = sum(sum(abs(Zy).^2));
            yx = sum(sum(abs(Zx).^2));
            yt = yx + yy;
            TE = yx/yt;
            TM = yy/yt;
            if TE > TM
                pol{ii} = 'TE';
            else
                pol{ii} = 'TM';
            end
        end

        results = cell(length(indices), 1);
        for ii  = 1:length(indices)
            temp = struct(...
                'wl', options(1).('wl'),...
                'neff', real(nProp(ii)),...
                'alpha', 4*pi/options(1).('wl')*imag(nProp(ii)),...
                'Polarization', pol{ii});
            results{ii} = temp;
        end

        % Cut off excess Modes
        if length(indices) > 2
            resultsTemp = cell(2, 1);
            for ii = 1:2
                maxval = 0;
               for ll = 1:length(results)
                   current = results{ll};
                   if current.('neff') > maxval
                       maxval = current.('neff');
                       I = ll;
                   end
               end
               resultsTemp{ii} = results{I};
               results(I) = [];
            end
            results = resultsTemp;
        end

        outerResults{iter} = results;
    end
end

