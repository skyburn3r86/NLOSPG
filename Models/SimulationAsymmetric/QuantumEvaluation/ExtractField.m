% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [Ep, Es, Ei] = ExtractField(model, varargin)
% Summary Here
%   Detailed Summary Here

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
        'hAir', {500, '[nm]'},...s
        'r33', {137, '[pm/V]'}, ...
        'lmin', {1100, '[nm]'}, ...
        'lmax', {2300, '[nm]'}, ...
        'Nl', {25, ' '}, ...
        'mat', {'', 'Si'}, ...
        'wl', {1550, '[nm]'}, ...
        'hRidge', {20, '[nm]'}, ...
        'dSlot', {100, '[nm]'}, ...
        'filename', {'', ''});

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
          warning('%s is not a recognized parameter name',inpName)
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
    if length(options(1).('wl')) > 1
        wl = options(1).('wl');
    else
        wl = linspace(options(1).('lmin'), options(1).('lmax'), options(1).('Nl'));
    end

    indices = zeros(1, 3);
    for ii = 1:length(wl)
        for kk = 1:length(wl)
            for ll = 1:length(wl)
                omegas = 2*pi*c_const/wl(ii);
                omegai = 2*pi*c_const/wl(kk);
                omegap = 2*pi*c_const/wl(ll);
                if  abs(omegap - omegai - omegas) < 2*pi*c_const/600e-6 && ii == kk
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
    
    modes = sortingModes(model); 
    neffp = modes{1}.TE(3);
    neffs = modes{2}.TE(1); 
    neffi = modes{2}.TE(1); 
    neffT = [neffs, neffi, neffp]; 
    % Extract fields as function of energy conservating modes
    for ii = 1:n
        for kk = 1:m
            iter = indices(ii, kk);
            temp = mpheval(model, 'ewfd.neff', 'dataset', 'dset2', 'outersolnum', iter);
            neff = temp.('d1');
            neff = neff(:, 1);
            [~, sI] = sort(real(neff));
            neff = neff(sI);
            
            temp = mpheval(model, 'ewfd.Ex', 'dataset', 'dset2', 'outersolnum', iter);
            Ex = temp.('d1');
            cordsx = temp.('p');

            temp = mpheval(model, 'ewfd.Ey', 'dataset', 'dset2', 'outersolnum', iter);
            Ey = temp.('d1');
            cordsy = temp.('p');

            temp = mpheval(model, 'ewfd.Ez', 'dataset', 'dset2', 'outersolnum', iter);
            Ez = temp.('d1');
            cordsz = temp.('p');

            temp = mpheval(model, 'ewfd.Hx', 'dataset', 'dset2', 'outersolnum', iter);
            Hx = temp.('d1');
            cordsvx = temp.('p');

            temp = mpheval(model, 'ewfd.Hy', 'dataset', 'dset2', 'outersolnum', iter);
            Hy = temp.('d1');
            cordsvy = temp.('p');

            temp = mpheval(model, 'ewfd.Hz', 'dataset', 'dset2', 'outersolnum', iter);
            Hz = temp.('d1');
            cordsvz = temp.('p');

            temp = mpheval(model, 'ewfd.nxx', 'dataset', 'dset2', 'outersolnum', iter);
            nRefr = temp.('d1');
            nRefr = nRefr(sI, :);

            Ex = Ex(sI, :);
            Ey = Ey(sI, :);
            Ez = Ez(sI, :);
            Hx = Hx(sI, :);
            Hy = Hy(sI, :);
            Hz = Hz(sI, :);

            x = cordsx(1, :);
            y = cordsx(2, :);

            n = 200;

            x_edge=linspace(min(x),max(x),n);
            y_edge=linspace(min(y),max(y),n);
            [X,Y]=meshgrid(x_edge,y_edge);
            sI = find(abs(neff - neffT(kk)) < 1e-4); 
            % TO DO: Only Consider important relevant modes.
            for zz = 1:length(sI)
                % Map Fields onto symmetric coordinate system.

                x = cordsx(1, :);
                y = cordsx(2, :);
                ux=griddata(x,y,Ex(sI(zz), :),X,Y);
                II = isnan(ux);
                ux(II) = 0;

                x = cordsy(1, :);
                y = cordsy(2, :);
                uy=griddata(x,y,Ey(sI(zz), :),X,Y);
                II = isnan(uy);
                uy(II) = 0;

                x = cordsz(1, :);
                y = cordsz(2, :);
                uz=griddata(x,y,Ez(sI(zz), :),X,Y);
                II = isnan(uz);
                uz(II) = 0;

                x = cordsvx(1, :);
                y = cordsvx(2, :);
                vx=griddata(x,y,Hx(sI(zz), :),X,Y);
                II = isnan(vx);
                vx(II) = 0;

                x = cordsvy(1, :);
                y = cordsvy(2, :);
                vy=griddata(x,y,Hy(sI(zz), :),X,Y);
                II = isnan(vy);
                vy(II) = 0;

                x = cordsvz(1, :);
                y = cordsvz(2, :);
                vz=griddata(x,y,Hz(sI(zz), :),X,Y);
                II = isnan(vz);
                vz(II) = 0;

                dx = X(1, 2) - X(1, 1);
                dy = Y(2, 1) - Y(1, 1);
                
                % Write Chi2 Matrix, for the objects
                ymin = options(1).('hSubstrate')+options(1).('hRidge');
                ymax = ymin + options(1).('hWG');
                xmin = (options(1).('wSim') - options(1).('dSlot'))/2;
                xmax = (options(1).('wSim') + options(1).('dSlot'))/2;
                
                x = x_edge(1,:);
                y = y_edge(1,:);
                chi2 = 0.1e-9;
                chi2_333 = zeros(length(x), length(y));
                for oo = 1:length(x)
                    for ll = 1:length(y)
                        if (x(oo) > xmin && x(oo) < xmax) && (y(ll) > ymin && y(ll) < ymax)
                            chi2_333(oo, ll) = 1;
                        end
                    end
                end
                ymin = options(1).('hSubstrate')+options(1).('hRidge')+options(1).('hWG');
                ymax = ymin + options(1).('hOrganic');
                xmin = 0; 
                xmax = options(1).('wSim');
                
                for oo = 1:length(x)
                    for ll = 1:length(y)
                        if (x(oo) > xmin && x(oo) < xmax) && (y(ll) > ymin && y(ll) < ymax)
                            chi2_333(oo, ll) = 1;
                        end
                    end
                end
                
                chi = chi2_333'*chi2;
                                    
                % Save Field Profile
                figure('Visible', 'off')
                surf(X*1e6,Y*1e6,abs(sqrt(conj(ux)*ux + conj(uy)*uy + conj(uz)*uz)))
                xlabel('x [um]')
                ylabel('y [um]')
                title([num2str(options(1).('wl')(iter)) 'nm: E for ' num2str(real(neff(ii))) ' + ' num2str(imag(neff(ii))) +'i'])
                saveas(gcf, ['./Figures/FieldProfiles/' options(1).('filename') '_' num2str(options(1).('wl')(iter)*1e9) 'nm_' num2str(real(neff(ii)), 4) '-i' num2str(abs(imag(neff(ii))), 4) '_Enorm.png'])

                % Asssign Fields
                if kk ~= 3
                    field = QuantumField(x_edge(1,:), y_edge(1,:), ux, uy, uz, vx, vy, vz, neff(sI(zz)), chi, wl(iter));
                else
                    field = ClassicalField(x_edge(1,:), y_edge(1,:), ux, uy, uz, vx, vy, vz, neff(sI(zz)), chi, wl(iter));
                end

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

    warning('on','all');
end