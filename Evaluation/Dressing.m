% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [x, y, Es, Ei, Ep] = Dressing(model, varargin)
%DRESSING Summary of this function goes here
%   Detailed explanation goes here
    % Parse Inputs
    warning('off','all');
    hbar = 1.0545718e-34;
    eps0 = 8.854e-12;
    mu0 = 4*pi*1e-7;
    c_const = 1/sqrt(eps0*mu0);

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
        'mat', {0, 'Si'}, ...
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

    % Find Energy Conservating Modes
    wl = linspace(options(1).('lmin'), options(1).('lmax'), options(1).('Nl'));
    indices = zeros(1, 3);
    for ii = 1:length(wl)
        for kk = 1:length(wl)
            for ll = 1:length(wl)
                omegas = 2*pi*c_const/wl(ii);
                omegai = 2*pi*c_const/wl(kk);
                omegap = 2*pi*c_const/wl(ll);
                if  abs(omegap - omegai - omegas) < 2*pi*c_const/300e-6
                   indices(end+1, :) = [ii, kk, ll];
                end
            end
        end
    end
    indices = indices(2:end, :);
    [n, m] = size(indices);

    Es = cell(length(n), 1);
    Ei = cell(length(n), 1);
    Ep = cell(length(n), 1);
    % Extract fields as function of energy conservating modes
    for ii = 1:n
        for kk = 1:m
            iter = indices(ii, kk);
            temp = mpheval(model, 'ewfd.neff', 'dataset', 'dset3', 'outersolnum', iter);
            neff = temp.('d1');
            neff = neff(:, 1);
            [~, sI] = sort(real(neff));
            neff = neff(sI);

            temp = mpheval(model, 'ewfd.Ex', 'dataset', 'dset3', 'outersolnum', iter);
            Ex = temp.('d1');
            cordsx = temp.('p');

            temp = mpheval(model, 'ewfd.Ey', 'dataset', 'dset3', 'outersolnum', iter);
            Ey = temp.('d1');
            cordsy = temp.('p');

            temp = mpheval(model, 'ewfd.Ez', 'dataset', 'dset3', 'outersolnum', iter);
            Ez = temp.('d1');
            cordsz = temp.('p');

            temp = mpheval(model, 'ewfd.nxx', 'dataset', 'dset3', 'outersolnum', iter);
            nRefr = temp.('d1');
            nRefr = nRefr(sI, :);

            Ex = Ex(sI, :);
            Ey = Ey(sI, :);
            Ez = Ey(sI, :);

            x = cordsx(1, :);
            y = cordsx(2, :);

            n = 100;

            x_edge=linspace(min(x),max(x),n);
            y_edge=linspace(min(y),max(y),n);
            [X,Y]=meshgrid(x_edge,y_edge);


            % TO DO: Only Consider important relevant modes.
            for zz = 1:length(sI)
                % Map Fields onto symmetric coordinate system.

                x = cordsx(1, :);
                y = cordsx(2, :);
                Zx=griddata(x,y,Ex(zz, :),X,Y);
                II = isnan(Zx);
                Zx(II) = 0;

                x = cordsy(1, :);
                y = cordsy(2, :);
                Zy=griddata(x,y,Ey(zz, :),X,Y);
                II = isnan(Zy);
                Zy(II) = 0;

                x = cordsy(1, :);
                y = cordsy(2, :);
                Zz=griddata(x,y,Ez(zz, :),X,Y);
                II = isnan(Zz);
                Zz(II) = 0;

                Zn = griddata(x, y, nRefr(ii, :), X, Y);
                II = isnan(Zn);
                Zn(II) = 0;
                dx = X(1, 2) - X(1, 1);
                dy = Y(2, 1) - Y(1, 1);

                Z = abs(sqrt(Zx.*conj(Zx) + Zy.*conj(Zy) + Zz.*conj(Zz)));
                mt = sum(sum(Z*dx*dy));
                mx = sum(sum(Z.*X*dx*dy));
                my = sum(sum(Z.*Y*dx*dy));
                % Select Relevant Modes
                xc = mx/mt;
                yc = my/mt;
                if ~((xc > (options(1).('wSim') - options(1).('wWG'))/2) && (xc < (options(1).('wSim') + options(1).('wWG'))/2) && ...
                    (yc > (options(1).('hSubstrate') - options(1).('hWG'))) && (yc < (options(1).('hSubstrate') + options(1).('hOrganic'))))
                    continue;
                end
                % Calculate Center of Energy
                % Select array elements to be evaluated
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
                pk = findpeaks(Z(0.9*length(Z), :));
                pk2 = findpeaks(Z(yi, xi));
                if ~((ratio > r && length(pk2) == 1) || (length(pk) == 1 && length(pk2) == 1))
                    continue
                end
                % Find Polarization
                yy = sum(sum(abs(Zy).^2));
                yx = sum(sum(abs(Zx).^2));
                yt = yx + yy;
                TE = yx/yt;
                TM = yy/yt;
                if TE > TM
                    disp('TE Mode. Ignoring')
                else
                    % Add fields into struct.
                    field = Field(Zx, Zy, Zz, 2*pi*c_const/wl(iter), neff(zz));
                    field = field.normalizeFields(x_edge, y_edge);
                    if kk == 1
                        Es{ii} = field;
                    elseif kk == 2
                        Ei{ii} = field;
                    elseif kk == 3
                        Ep{ii} = field;
                    end
                end

            end
        end
    end
    x = x_edge;
    y = y_edge;
    warning('on','all');
end

