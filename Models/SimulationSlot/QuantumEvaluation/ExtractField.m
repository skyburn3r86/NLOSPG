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

    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'OuterSolNums'
                OuterSolNums = varargin{ii+1};
            case 'SolNums'
                SolNums = varargin{ii+1};
        end
    end
    

    % Extract fields as function of energy conservating modes
    for ii = 1:length(OuterSolNums)
        iter = OuterSolNums(ii);
        neff = mphglobal(model, 'ewfd.neff', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums(ii));

        temp = mpheval(model, 'ewfd.Ex', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums(ii));
        Ex = temp.('d1');
        cordsx = temp.('p');

        temp = mpheval(model, 'ewfd.Ey', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums(ii));
        Ey = temp.('d1');
        cordsy = temp.('p');

        temp = mpheval(model, 'ewfd.Ez', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums(ii));
        Ez = temp.('d1');
        cordsz = temp.('p');

        temp = mpheval(model, 'ewfd.Hx', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums(ii));
        Hx = temp.('d1');
        cordsvx = temp.('p');

        temp = mpheval(model, 'ewfd.Hy', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums(ii));
        Hy = temp.('d1');
        cordsvy = temp.('p');

        temp = mpheval(model, 'ewfd.Hz', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums(ii));
        Hz = temp.('d1');
        cordsvz = temp.('p');

        x = cordsx(1, :);
        y = cordsx(2, :);

        n = 200;

        x_edge=linspace(min(x),max(x),n);
        y_edge=linspace(min(y),max(y),n);
        [X,Y]=meshgrid(x_edge,y_edge);

        % Map Fields onto symmetric coordinate system.

        x = cordsx(1, :);
        y = cordsx(2, :);
        ux=griddata(x,y,Ex,X,Y);
        II = isnan(ux);
        ux(II) = 0;

        x = cordsy(1, :);
        y = cordsy(2, :);
        uy=griddata(x,y,Ey,X,Y);
        II = isnan(uy);
        uy(II) = 0;

        x = cordsz(1, :);
        y = cordsz(2, :);
        uz=griddata(x,y,Ez,X,Y);
        II = isnan(uz);
        uz(II) = 0;

        x = cordsvx(1, :);
        y = cordsvx(2, :);
        vx=griddata(x,y,Hx,X,Y);
        II = isnan(vx);
        vx(II) = 0;

        x = cordsvy(1, :);
        y = cordsvy(2, :);
        vy=griddata(x,y,Hy,X,Y);
        II = isnan(vy);
        vy(II) = 0;

        x = cordsvz(1, :);
        y = cordsvz(2, :);
        vz=griddata(x,y,Hz,X,Y);
        II = isnan(vz);
        vz(II) = 0;

        dx = X(1, 2) - X(1, 1);
        dy = Y(2, 1) - Y(1, 1);

        % Parse Model Parameters
        hRidge = char(model.param.get('hRidge'));
        if contains(hRidge, '[um]')
            hRidge = str2double(erase(hRidge, ' [um]'))*1e-6;
        elseif contains(hRidge, '[nm]')
            hRidge = str2double(erase(hRidge, ' [nm]'))*1e-9;
        elseif contains(hRidge, '[m]')
            hRidge = str2double(erase(hRidge, ' [m]'));
        else
            error('hRidge: Unit not recognized')
        end

        hSubstrate = char(model.param.get('hSubstrate'));
        if contains(hSubstrate, '[um]')
            hSubstrate = str2double(erase(hSubstrate, ' [um]'))*1e-6;
        elseif contains(hSubstrate, '[nm]')
            hSubstrate = str2double(erase(hSubstrate, ' [nm]'))*1e-9;
        elseif contains(hSubstrate, '[m]')
            hSubstrate = str2double(erase(hSubstrate, ' [m]'));
        else
            error('hSubstrate: Unit not recognized')
        end

        hWG = char(model.param.get('hWG'));
        if contains(hWG, '[um]')
            hWG = str2double(erase(hWG, ' [um]'))*1e-6;
        elseif contains(hWG, '[nm]')
            hWG = str2double(erase(hWG, ' [nm]'))*1e-9;
        elseif contains(hWG, '[m]')
            hWG = str2double(erase(hWG, ' [m]'));
        else
            error('hWG: Unit not recognized')
        end

        hOrganic = char(model.param.get('hOrganic'));
        if contains(hOrganic, '[um]')
            hOrganic = str2double(erase(hOrganic, ' [um]'))*1e-6;
        elseif contains(hOrganic, '[nm]')
            hOrganic = str2double(erase(hOrganic, ' [nm]'))*1e-9;
        elseif contains(hOrganic, '[m]')
            hOrganic = str2double(erase(hOrganic, ' [m]'));
        else
            error('hOrganic: Unit not recognized')
        end

        dSlot = char(model.param.get('dSlot'));
        if contains(dSlot, '[um]')
            dSlot = str2double(erase(dSlot, ' [um]'))*1e-6;
        elseif contains(dSlot, '[nm]')
            dSlot = str2double(erase(dSlot, ' [nm]'))*1e-9;
        elseif contains(dSlot, '[m]')
            dSlot = str2double(erase(dSlot, ' [m]'));
        else
            error('dSlot: Unit not recognized')
        end

        wSim = char(model.param.get('wSim'));
        if contains(wSim, '[um]')
            wSim = str2double(erase(wSim, ' [um]'))*1e-6;
        elseif contains(wSim, '[nm]')
            wSim = str2double(erase(wSim, ' [nm]'))*1e-9;
        elseif contains(wSim, '[m]')
            wSim = str2double(erase(wSim, ' [m]'));
        else
            error('wSim: Unit not recognized')
        end

        % Write Chi2 Matrix, for the objects
        ymin = hSubstrate + hRidge; 
        ymax = ymin + hWG;
        xmin = (wSim - dSlot)/2;
        xmax = xmin + dSlot; 

        x = x_edge(1,:);
        y = y_edge(1,:);
        chi2 = 0.1e-9;
        chi2_333 = zeros(length(x), length(y));
        for oo = 1:length(x)
            for ll = 1:length(y)
                if (x(oo) > xmin && x(oo) < xmax) && (y(ll) > ymin && y(ll) < ymax)
                    chi2_333(ll, oo) = 1;
                end
            end
        end

        ymin = ymax;
        ymax = ymin + hOrganic;
        xmin = 0; 
        xmax = wSim;

        for oo = 1:length(x)
            for ll = 1:length(y)
                if (x(oo) > xmin && x(oo) < xmax) && (y(ll) > ymin && y(ll) < ymax)
                    chi2_333(ll, oo) = 1;
                end
            end
        end

        chi = chi2_333*chi2;
        
        % extracting group refractive index via Sum(Energy)/Sum(PowerFlow) -
        % All-plasmonic Mach-Zehnder Modulator Haffner et al. Nature Photonics
        % (2015)
        [mode_energy.value, mode_energy.unit] = mphint2(model,...
            '(eps0*(ewfd.nxx^2*abs(ewfd.Ex)^2+ewfd.nyy^2*abs(ewfd.Ey)^2+ewfd.nzz^2*abs(ewfd.Ez)^2))',...
            'surface', 'dataset', 'dset2', 'solnum', SolNums(ii), 'outersolnum', iter);   
        [mode_power_flow.value, mode_power_flow.unit] = mphint2(model,...
            '(ewfd.Ex*conj(ewfd.Hy)-conj(ewfd.Ey)*(ewfd.Hx))',...
            'surface', 'dataset', 'dset2', 'solnum', SolNums(ii), 'outersolnum', iter); 
        ng = real(3e8*mode_energy.value/mode_power_flow.value); 

        % Asssign Fields
        if ii ~= 1
            field = QuantumField(x_edge(1,:), y_edge(1,:), ux, uy, uz, vx, vy, vz, neff, real(ng), chi, 2620e-9);
        else
            field = ClassicalField(x_edge(1,:), y_edge(1,:), ux, uy, uz, vx, vy, vz, neff, real(ng), chi, 1310e-9);
        end

        if ii == 1
            Ep = field;
        elseif ii == 2
            Ei = field;
            Es = field;
        end

    end

    warning('on','all');
end