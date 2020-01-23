% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [P, P2] = SPDC_QS(x, y, chi, Es, Ei, Ep, omegas, neffs, omegai, neffi, omegap, neffp, l, v, Ppump)
%SPDC_QS Summary of this function goes here
%   Detailed explanation goes here
    % after Fiorentino, M., Spillane, S. M., Beausoleil, R. G., Roberts, T. D., Battle, P., & Munro, M. W. (2007).
    % Spontaneous parametric down-conversion in periodically poled KTP waveguides and bulk crystals. Optics Express, 15(12), 7479. https://doi.org/10.1364/oe.15.007479
    % Physical Constants
    phase = "QuasiPhase";
    losses = false;
    eps0 = 8.854e-12;           % F/m
    mu0 = 4*pi*1e-7;            % H/m
    c_const = 1/sqrt(mu0*eps0); % m/s
    hbar = 1.10545718e-34;      % Joules * s


    dx = abs(x(2) - x(1));
    dy = abs(y(2) - y(1));

    z = linspace(0, l, 3000);
    dz = abs(z(2)-z(1));
    g = sqrt(2/(eps0*neffs^2*neffi^2*neffp))*sqrt(Ppump)/l*sqrt(omegas*omegai);
    It = l/v;
    Eys = Es.Ey;
    Eyi = Ei.Ey;
    Eyp = Ep.Ey;
    if strcmp(phase, "Perfect")
        dk = 0;
    else
        dk = real(neffp*omegap/c_const - neffs*omegas/c_const - neffi*omegai/c_const);
        lc = pi/dk;
        T = 2*lc;
    end

    Iap = 1/sqrt(l)^3*dx*dy*sum(sum(chi.*conj(Eys).*conj(Eyi).*Eyp));
    Iam = 1/sqrt(l)^3*dx*dy*sum(sum((-chi).*conj(Eys).*conj(Eyi).*Eyp));

    % Preperation of Hamiltonian
    NPhotons = 30;
    n = 1:NPhotons;

    % Operator Preparation
    a = diag(sqrt(n), 1);
    ad = diag(sqrt(n), 1);
    ad = ad';
    % Initialize State Vector
    v1 = zeros(NPhotons+1, 1);
    v1(1) = 1;
    vF = zeros(NPhotons+1, 1);
    vF2 = vF;
    vF2(3) = 1;
    vF(2) = 1;
    if losses == false
        neffs = real(neffs);
    end
    Il = zeros(length(v1), length(z)) + 1i*zeros(length(v1), length(z));
    Il(:, 1) = v1;
    if dk ~= 0
        if strcmp(phase, "QuasiPhase")
            for ii = 2:length(z)
               if mod(z(ii), T) <= lc
                   Il(:, ii-1) = Il(:, ii-1)/sqrt(Il(:, ii-1)'*Il(:, ii-1));
                   Il(:, ii) = Il(:, ii-1) + Iap*exp(1i*dk*z(ii))*dz*ad*Il(:, ii-1)*exp(-imag(neffs)*omegas/c_const*dz);
               elseif mod(z(ii), T) > lc
                   Il(:, ii-1) = Il(:, ii-1)/sqrt(Il(:, ii-1)'*Il(:, ii-1));
                   Il(:, ii) = Il(:, ii-1) + Iam*exp(1i*dk*z(ii))*dz*ad*Il(:, ii-1)*exp(-imag(neffs)*omegas/c_const*dz);
               end
            end
        elseif strcmp(phase, "None")
            for ii = 2:length(z)
                Il(:, ii) = Il(:, ii-1) + Iap*exp(1i*dk*z(ii))*dz*ad*Il(:, ii-1)/sqrt(Il(:, ii-1)'*Il(:, ii-1));
            end
        end
    else
        for ii = 2:length(z)
            Il(:, ii) = Il(:, ii-1) + Iap*dz*ad*Il(:, ii-1)/sqrt(Il(:, ii-1)'*Il(:, ii-1));
        end
    end
    P = -1i*g*It*vF'*Il*hbar;
    P2 = -1i*g*It*vF2'*Il*hbar;
end

