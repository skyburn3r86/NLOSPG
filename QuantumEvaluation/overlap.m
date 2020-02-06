% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function g0 = overlap(Es, Ep, Ei)
    x = Es.x;
    y = Ei.y;
    chi = Es.chi2r;
%     int = chi*conj(Es.uy)*conj(Ei.uy)*Ep.uy;
    int = (Es.uy)*conj(Ep.uy)*200e-6;
    A = sqrt(Es.omega*Ei.omega)/(2*real(Es.neff)*real(Ei.neff))*sqrt(Ep.hbar*Ep.omega/(2*real(Ep.neff)^2*Ep.eps0));       
    overlap = trapz(y, trapz(x, int, 2));
    g0 = A*overlap;
end