% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function g0 = overlap(Es, Ep, Ei)
    x = Es.x;
    y = Ei.y;
    chi = Es.chi2r;
    int = chi*conj(Es.uy)*conj(Ei.uy)*Ep.uy;
    A = Es.hbar*sqrt(Es.omega*Ei.omega)/(2*Es.eps0*real(Es.neff)*real(Ei.neff))*sqrt(Ep.hbar*Ep.omega/(2*real(Ep.neff)^2*Ep.eps0));       % Missing Factor 1/sqrt(L) due to field normalization
    overlap = trapz(y, trapz(x, int, 2));
    g0 = A*overlap;
end