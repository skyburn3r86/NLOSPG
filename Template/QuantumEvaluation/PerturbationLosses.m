% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [field] = PerturbationLosses(field, t)
    H = field.a;
    hbar = field.hbar;
    l = (1i*field.chi1i*field.hbar*field.omega);
    
    NPhotons = field.NPhotons + 1;
    gamma0 = zeros(NPhotons, 1);
    gamma1 = zeros(NPhotons, 1);
    gamma2 = zeros(NPhotons, 1);
    gamma3 = zeros(NPhotons, 1);
    gamma4 = zeros(NPhotons, 1); 
    % Find Initial Conditions
    for kk = 1:NPhotons
        psik = zeros(NPhotons, 1);
        psik(kk) = 1;
        gamma0(kk) = psik'*full(field.psi);
    end
    for kk = 1:NPhotons
        psik = zeros(NPhotons, 1);
        psik(kk) = 1;
        gamma1(kk) = -1i/hbar*(psik'*(H*full(field.psi)))*t;                                                        % First Perturbation Coefficient
        for ii = 1:NPhotons
            psin = zeros(NPhotons, 1);
            psin(ii) = 1;
            gamma2(kk) = gamma2(kk)-1/(2*hbar^2)*(psin'*H*full(field.psi))*(psik'*H*psin)*t^2;                      % Second Perturbation Coefficient
            for ll = 1:NPhotons
               psil = zeros(NPhotons, 1); 
               psil(ll) = 1; 
               gamma3(kk) = gamma3(kk) -1i/(6*hbar^3)*(psil'*H*full(field.psi))*(psin'*H*psil)*(psik'*H*psin)*t^3;  % Third Perturbation Coefficient
               for nn = 1:NPhotons
                  gamma4(kk) = gamma4(kk);
               end
            end
        end
    end
    gamma = gamma0 + l*gamma1 + l^2*gamma2 + l^3*gamma3;
    field.psi = gamma; 
    field.psi = 1/sqrt(field.psi'*field.psi)*field.psi; 
end