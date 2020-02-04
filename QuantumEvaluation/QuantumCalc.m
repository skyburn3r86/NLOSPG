% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [P, Es, Ei, Ep] = QuantumCalc(Es, Ep, Ei, z)

    [Es, Ei, Ep] = Propagate(Es, Ei, Ep, z); 
    P = zeros(length(z), 1); 
    for ii = 1:length(z)
       P(ii) = power(Ep{ii});
    end

end

% function Ep = propagate(Ep, z)
%     x = Ep.x;
%     y = Ep.y;
%     chi1i = Ep.chi1i;
%     omega = Ep.omega;
%     eps0 = Ep.eps0;
%     epsr = real(Ep.neff)^2;
%     
%     if length(z) > 1
%         dz = z(2) - z(1);
%     else
%         dz = 0;
%     end
%     
%     A = epsr/(Ep.c_const)*chi1i*omega;                     % Right solution would be: Epsr*c_const/n_G
%     sigma = sparse(diag(ones(Ep.NPhotons, 1), 1)); 
%     
%     for ii = 1:length(z)
%         N = full(Ep.psi'*Ep.ad*Ep.a*Ep.psi);
%         Ep.psi = expA(A*dz*Ep.a)*Ep.psi; 
%         B = sqrt(Ep.psi'*Ep.psi); 
%         Ep.psi = 1/B*Ep.psi; 
%     end
% end
% 

function P = power(field)
    dx = field.x(2) - field.x(1); 
    dy = field.y(2) - field.y(1); 
    N = field.psi'*field.ad*field.a*field.psi; 
    S = sqrt(field.eps0/field.mu0); 
    Sz = real(field.hbar*field.omega/(2*field.eps0*sqrt(real(field.neff)))*(field.ux*conj(field.vy) - field.uy*conj(field.vx))*(N));  
    P = S*dx*dy*sum(sum(Sz));
end