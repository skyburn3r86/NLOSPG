% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [Es, Ei, Ep, g0] = Propagate(Es, Ei, Ep, z)
    % Specification of the Type of Phase Matching
    phase = "Perfect";
    
    Estemp = cell(length(z), 1); 
    Eitemp = cell(length(z), 1); 
    Eptemp = cell(length(z), 1); 
    
    for ii = 1:length(z)
       Estemp{ii} = Es; 
       Eitemp{ii} = Ei;
       Eptemp{ii} = Ep;
    end
    
    % Constants
    eps0 = 8.854e-12;           % F/m
    mu0 = 4*pi*1e-7;            % H/m
    c_const = 1/sqrt(eps0*mu0); % m/s
    hbar = 1.10545718e-34;      % Js

    % Get Coordinate system
    x = Ep.x;
    y = Ep.y;

    % Momentum Mismatch
    dk = real(Ep.neff)*Ep.omega/c_const - (real(Es.neff)*Es.omega/c_const + real(Ei.neff)*Ei.omega/c_const);
    lc = pi/dk;
    T = 2*lc;

    if length(z) > 1
        dz = z(2) - z(1);
    else
        dz = 0;
    end
    % Loss rate [1/m]
    gLp = Ep.chi1i/real(Ep.neff)*Ep.omega/c_const;
    g0 = overlap(Es, Ep, Ei);

    for ii = 2:length(z)
        N = Eptemp{ii-1}.N;
        % Non-Linear Process
        if strcmp(phase, "QuasiPhase")          % Periodically switch poling direction
            if mod(z(ii), T) <= lc
                Estemp{ii}.psi = Estemp{ii-1}.psi + g0*exp(1i*dk*z(ii))*dz/c_const*Es.ad*Estemp{ii-1}.psi*sqrt(N);  % Should be dz
                Eitemp{ii}.psi = Eitemp{ii-1}.psi + g0*exp(1i*dk*z(ii))*dz/c_const*Ei.ad*Eitemp{ii-1}.psi*sqrt(N);  % 
            else
                Estemp{ii}.psi = Estemp{ii-1}.psi - g0*exp(1i*dk*z(ii))*dz/c_const*Es.ad*Estemp{ii-1}.psi*sqrt(N);
                Eitemp{ii}.psi = Eitemp{ii-1}.psi - g0*exp(1i*dk*z(ii))*dz/c_const*Ei.ad*Eitemp{ii-1}.psi*sqrt(N);
            end

        elseif strcmp(phase, "Perfect")         % No Phase Mismatch
            Estemp{ii}.psi = Estemp{ii-1}.psi + g0*dz/c_const*Estemp{ii-1}.ad*Estemp{ii-1}.psi*sqrt(N);
            Eitemp{ii}.psi = Eitemp{ii-1}.psi + g0*dz/c_const*Eitemp{ii-1}.ad*Eitemp{ii-1}.psi*sqrt(N);

        elseif strcmp(phase, "None")            % Phase Mismatch non-compensated
            Estemp{ii}.psi = Estemp{ii-1}.psi + g0*exp(1i*dk*dz)*dz/c_const*Es.ad*Estemp{ii-1}.psi*sqrt(N);
            Eitemp{ii}.psi = Eitemp{ii-1}.psi + g0*exp(1i*dk*dz)*dz/c_const*Ei.ad*Eitemp{ii-1}.psi*sqrt(N);
        end

        Estemp{ii}.psi = 1/sqrt(Estemp{ii}.psi'*Estemp{ii}.psi)*Estemp{ii}.psi;         % Normalization
        Eitemp{ii}.psi = 1/sqrt(Eitemp{ii}.psi'*Eitemp{ii}.psi)*Eitemp{ii}.psi;         % Normalization

        % Losses
        % Losses Pump
        Eptemp{ii} = Eptemp{ii}.setPower(Eptemp{ii-1}.Pin*exp(gLp*dz));                 % Propagation

        % Losses Signal
        Estemp{ii} = losses(Estemp{ii}, dz, 'type', 'Beamsplitter');                     % Propagation

        % Losses Idler
        Eitemp{ii} = losses(Eitemp{ii}, dz, 'type', 'Beamsplitter');                     % Propagation

    end
    Es = Estemp;
    Ep = Eptemp; 
    Ei = Eitemp; 
end


function s = expA(x)
    s = zeros(size(x)); 
    N = length(x)/2;
    for ii = 0:N
        s = s + 1/factorial(ii)*x^ii;
    end
end