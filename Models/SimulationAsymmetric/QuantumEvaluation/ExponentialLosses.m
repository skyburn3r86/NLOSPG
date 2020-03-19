% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [field] = ExponentialLosses(field, z)
% Based on Rana, F. (2009). Chapter 9: Loss in Quantum Optics. Compression of operator
%   Detailed Summary Here
    gamma = -1/(2*field.neff)*field.chi1i*field.omega/field.c_const;
    field.a = field.a*exp(-2*gamma*z);
    field.ad = field.ad*exp(-2*gamma*z);
end