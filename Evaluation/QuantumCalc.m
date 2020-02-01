% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [g0, P, Ex, Ey, Ez, Hx, Hy, Hz] = QuantumCalc(Es, Ep, Ei, z)
% Summary Here
%   Detailed Summary Here
    %g0 = overlap(Es, Ep, Ei); 
    g0 = 0; 
    Ep = propagate(Ep, z); 
    [Ex, Ey, Ez, Hx, Hy, Hz] = Eexp(Ep); 
    P = power(Ex, Ey, Ez, Hx, Hy, Hz, Ep); 
end

function g0 = overlap(Es, Ep, Ei)
    x = Es.x;
    y = Ei.y;
    chi = Es.chi2r;
    int = chi*conj(Es.uy)*conj(Ei.uy)*Ep.uy;
    A = Es.hbar*sqrt(Es.omega*Ei.omega)/(2*Es.eps0*real(Es.neff)*real(Ei.neff))*sqrt(Ep.Pp/(Ep.c_const*real(Ep.neff)*Ep.eps0));
    overlap = trapz(y, trapz(x, int, 2));
    g0 = A*overlap;
end

function Ep = propagate(Ep, z)
    x = Ep.x;
    y = Ep.y;
    chi1i = Ep.chi1i;
    omega = Ep.omega;
    eps0 = Ep.eps0;
    epsr = real(Ep.neff)^2;
    if length(z) > 1
        dz = z(2) - z(1);
    else
        dz = 0;
    end
    A = sqrt(epsr)*chi1i*omega/(Ep.c_const/real(Ep.neff));
    sigma = sparse(diag(ones(Ep.NPhotons, 1), 1)); 
    for ii = 1:length(z)
%         Ep.psi = (eye(size(sigma))+A*(ii*dz)*Ep.a + 0.5*(A*(ii*dz))^2*Ep.a*Ep.a)*Ep.psi; 
%         Ep.psi = (eye(size(sigma))+A*dz*Ep.a + 0.5*(A*dz)^2*Ep.a*Ep.a + 1/6*(A*dz)^3*Ep.a*Ep.a*Ep.a)*Ep.psi; 
        Ep.psi = (eye(size(sigma))+A*dz*sigma)*Ep.psi; 
        Ep.psi = Ep.psi/sqrt((Ep.psi'*Ep.ad*Ep.a*Ep.psi));
    end
end
% Moving back to the classical Domain:
% Expectation value of the Electric Field 
function [Ex, Ey, Ez, Hx, Hy, Hz] = Eexp(field)
    psit = ones(length(field.psi), 1); 
    
%     Ex = 1i*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.ux*(psit'*field.a*field.psi) - conj(field.ux)*(psit'*field.ad*field.psi));
%     Hx = 1i/field.omega*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.vx*(psit'*field.a*field.psi) - conj(field.vx)*(psit'*field.ad*field.psi));
%     Ey = 1i*field.uy*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.uy*(psit'*field.a*field.psi) - conj(field.uy)*(psit'*field.ad*field.psi));
%     Hy = 1i/field.omega*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.vy*(psit'*field.a*field.psi) - conj(field.vy)*(psit'*field.ad*field.psi));
%     Ez = 1i*field.uz*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.uz*(psit'*field.a*field.psi) - conj(field.uz)*(psit'*field.ad*field.psi));
%     Hz = 1i/field.omega*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.vz*(psit'*field.a*field.psi) - conj(field.vz)*(psit'*field.ad*field.psi));
    
    
    Ex = 1i*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.ux*(psit'*field.a*field.psi));
    Hx = 1i/field.omega*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.vx*(psit'*field.a*field.psi));
    Ey = 1i*field.uy*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.uy*(psit'*field.a*field.psi));
    Hy = 1i/field.omega*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.vy*(psit'*field.a*field.psi));
    Ez = 1i*field.uz*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.uz*(psit'*field.a*field.psi));
    Hz = 1i/field.omega*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*(field.vz*(psit'*field.a*field.psi));
end

function P = power(Ex, Ey, Ez, Hx, Hy, Hz, field)
    dx = field.x(2) - field.x(1); 
    dy = field.y(2) - field.y(1); 
    Sz = real(Ex*conj(Hy) - Ey*conj(Hx)); 
    P = dx*dy*abs(sum(sum(Sz)));
end