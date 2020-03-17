% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [field] = ScatteringLosses(field, z)
% Scattering into a lower order mode with a certain probability
%   Detailed Summary Here
    gamma = 0.1;
    psi1 = field.a*field.psi;
    psi1 = psi1/sqrt((psi1'*psi1));
    field.psi = field.psi*(1-gamma) + psi1*gamma;
    field.psi = field.psi/sqrt(field.psi'*field.psi);

end