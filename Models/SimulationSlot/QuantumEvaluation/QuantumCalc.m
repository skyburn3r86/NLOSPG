% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [P, Es, Ei, Ep, g0] = QuantumCalc(Es, Ep, Ei, z)

    [Es, Ei, Ep, g0] = Propagate(Es, Ei, Ep, z);
    P = zeros(length(z), 1); 
    for ii = 1:length(z)
       P(ii) = power(Ep{ii});
    end
end

function P = power(field)
    if isa(field, 'ClassicalField')
        P = field.Pin;
    else
        dx = field.x(2) - field.x(1);
        dy = field.y(2) - field.y(1);
        N = field.psi'*field.ad*field.a*field.psi;
        S = sqrt(field.eps0/field.mu0);
        Sz = real(field.hbar*field.omega/(2*field.eps0*sqrt(real(field.neff)))*(field.ux*conj(field.vy) - field.uy*conj(field.vx))*(N));
        P = S*dx*dy*sum(sum(Sz));
    end
end