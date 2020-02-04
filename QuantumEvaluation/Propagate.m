% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [Es, Ei, Ep] = Propagate(Es, Ei, Ep, z)
    % Specification of the Type of Phase Matching
    debug = true; 
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
    gLp = real(Ep.neff)^2/(c_const)*Ep.chi1i*Ep.omega;                     % Right solution would be: Epsr*c_const/n_G
    gLi = real(Ei.neff)^2/(c_const)*Ei.chi1i*Ei.omega;
    gLs = real(Es.neff)^2/(c_const)*Es.chi1i*Es.omega;

    g0 = overlap(Es, Ep, Ei);
    N = 1e24;
    for ii = 2:length(z)
        % Non-Linear Process
        if strcmp(phase, "QuasiPhase")          % Periodically switch poling direction
            if mod(z(ii), T) <= lc
                Estemp{ii}.psi = Estemp{ii-1}.psi + g0*exp(1i*dk*dz)*dz*Es.ad*Estemp{ii-1}.psi*sqrt(N);
                Eitemp{ii}.psi = Eitemp{ii-1}.psi + g0*exp(1i*dk*dz)*dz*Ei.ad*Eitemp{ii-1}.psi*sqrt(N);
            else
                Estemp{ii}.psi = Estemp{ii-1}.psi - g0*exp(1i*dk*dz)*dz*Es.ad*Estemp{ii-1}.psi*sqrt(N);
                Eitemp{ii}.psi = Eitemp{ii-1}.psi - g0*exp(1i*dk*dz)*dz*Ei.ad*Eitemp{ii-1}.psi*sqrt(N);
            end

        elseif strcmp(phase, "Perfect")         % No Phase Mismatch
            Estemp{ii}.psi = Estemp{ii-1}.psi + g0*dz*Es.ad*Estemp{ii-1}.psi*sqrt(N);
            Eitemp{ii}.psi = Eitemp{ii-1}.psi + g0*dz*Ei.ad*Eitemp{ii-1}.psi*sqrt(N);

        elseif strcmp(phase, "None")            % Phase Mismatch non-compensated
            Estemp{ii}.psi = Estemp{ii-1}.psi + g0*exp(1i*dk*dz)*dz/c_const*Es.ad*Estemp{ii-1}.psi*sqrt(N);
            Eitemp{ii}.psi = Eitemp{ii-1}.psi + g0*exp(1i*dk*dz)*dz/c_const*Ei.ad*Eitemp{ii-1}.psi*sqrt(N);
        end

        Estemp{ii}.psi = 1/sqrt(Estemp{ii}.psi'*Estemp{ii}.psi)*Estemp{ii}.psi;         % Normalization
        Eitemp{ii}.psi = 1/sqrt(Eitemp{ii}.psi'*Eitemp{ii}.psi)*Eitemp{ii}.psi;         % Normalization

        % Losses
        % Losses Pump
        Eptemp{ii}.psi = expA(dz*gLp*Ep.a)*Eptemp{ii-1}.psi;              % Propagation
        Eptemp{ii}.psi = 1/sqrt(Eptemp{ii}.psi'*Eptemp{ii}.psi)*Eptemp{ii}.psi;         % Normalization

        % Losses Signal
        Estemp{ii}.psi = expA(dz*gLs*Es.a)*Estemp{ii}.psi;                            % Propagation
        Estemp{ii}.psi = 1/sqrt(Estemp{ii}.psi'*Estemp{ii}.psi)*Estemp{ii}.psi;         % Normalization

        % Losses Idler
        Eitemp{ii}.psi = expA(dz*gLi*Ei.a)*Eitemp{ii}.psi;                            % Propagation
        Eitemp{ii}.psi = 1/sqrt(Eitemp{ii}.psi'*Eitemp{ii}.psi)*Eitemp{ii}.psi;         % Normalization

    end
    Es = Estemp;
    Ep = Eptemp; 
    Ei = Eitemp; 
end

function g0 = overlap(Es, Ep, Ei)
    x = Es.x;
    y = Ei.y;
    chi = Es.chi2r;
    int = chi*conj(Es.uy)*conj(Ei.uy)*Ep.uy;
    A = Es.hbar*sqrt(Es.omega*Ei.omega)/(2*Es.eps0*real(Es.neff)*real(Ei.neff))*sqrt(Ep.hbar*Ep.omega/(2*real(Ep.neff)^2*Ep.eps0));       % Missing Factor 1/sqrt(L) due to field normalization
    overlap = trapz(y, trapz(x, int, 2));
    g0 = A*overlap;
end

function s = expA(x)
    s = zeros(size(x)); 
    N = length(x)/2;
    for ii = 0:N
        s = s + 1/factorial(ii)*x^ii;
    end
end