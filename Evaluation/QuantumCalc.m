% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [x, y] = QuantumCalc(Es, Ep, Ei)
% Summary Here
%   Detailed Summary Here


end

function g0 = overlap(Es, Ep, Ei)
    x = Es.x;
    y = Ey.y;
    chi = Es.chi;
    int = chi*conj(Es.uy)*conj(Ei.uy)*Ep.uy;
    A = Es.hbar*sqrt(Es.omega*Ei.omega)/(2*Es.eps0*real(Es.neff)*real(Ei*neff))*sqrt(Ep.Pp/(Ep.c_const*real(Ep.neff)*Ep.eps0))
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
    A = -eps0/epsr*chi1i*omega*z(end)/(Ep.c_const/real(Ep.neff));
    Ep.psi = Ep.psi + A*Ep.psi
end

function Ex, Ey, Ez = Eexp(field)
    Ex = i*field.ux*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*field.psi'*ad*a*field.psi;
    Ey = i*field.uy*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*field.psi'*ad*a*field.psi;
    Ey = i*field.uz*sqrt(field.hbar*field.omega/(2*sqrt(real(field.neff))*field.eps0))*field.psi'*ad*a*field.psi;
end

function P = power(Ex, Ey, Ez)
    P = 1/377*(Ex*conj(Ex) + Ey*conj(Ey) + Ez*conj(Ez))
end