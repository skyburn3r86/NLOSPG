% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function field = Beamsplitter(field, dz, zmax)
% Summary Here
%   Detailed Summary Here
    eta = -1/(2*field.neff)*field.chi1i;
    rho = field.psi*field.psi';
    coeffs = zeros(field.NPhotons + 1, 1);
    results = zeros(field.NPhotons + 1, field.NPhotons + 1); 
    for ii = 1:min(max(find(field.psi > 0)) + 1, field.NPhotons + 1)               % Can be accelerated by considering only the photons up to occupation and not beyond
        psik = zeros(field.NPhotons + 1, 1);
        psik(ii) = 1;
        coeffs(ii) = psik'*field.psi;                                       % Fish out the Probability coefficients coefficients of being in state k: gamma = <k|psi>
        p = zeros(field.NPhotons + 1, 1);
        for kk = 0:ii-1
            psil = zeros(field.NPhotons + 1, 1);                            
            psil(ii - kk) = 1;                                              % Trace out ii - kk photons, i.e. reduce the occupation number of the highest state by kk.
            p = p + nchoosek(ii, kk)*eta^kk*(1-eta)^(ii-kk)*psil;           % Calculate the probability of said transition
        end
        results(:, ii) = p;  
    end
    psi = results*coeffs;                                                   % Obtain only the relevant parts of the solution. All artificial solutions will be removed
    psi = psi/sqrt(psi'*psi);                                               % Renormalize
    field.psi = psi;                                                        % Write
end