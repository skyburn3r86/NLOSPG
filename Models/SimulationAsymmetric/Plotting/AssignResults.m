% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [SimRes] = AssignResults(results, hWG, hOrganic)
%ASSIGNRESULTS Summary of this function goes here
%   Detailed explanation goes here

    % Collect Data from results
    c_const = 299792458;
    save = true;
    close all
    flat = vertcat(results{:});

    wl = [];
    polarization = cell(length(flat),1);
    nA = [];
    kA = [];
    for ii = 1:length(flat)
        currentelement = flat{ii};
        lambda = currentelement.('wl');
        pol = currentelement.('Polarization');
        n = currentelement.('neff');
        alpha = currentelement.('alpha');
        wl = [wl lambda];
        polarization{ii} = pol;
        nA = [nA n];
        kA = [kA alpha];
    end
    [wl, I] = sort(wl);
    polarization = polarization(I);
    nA = nA(I);
    kA = kA(I);
    indices = [];
    ii = 1;
    while ii <= length(wl)-1
       if wl(ii) == wl(ii+1)
           ii = ii + 2;
           continue
       else
          indices = [indices ii];
          ii = ii + 1;
          if ii == length(wl)
              indices = [indices ii];
          end
       end
    end


    for ii = 1:length(indices)
        wl = [wl(1:indices(ii)) wl(indices(ii):end)];
        if strcmp(polarization(indices(ii)+(ii-1)), 'TM')
            polarization = {polarization{1:indices(ii)+(ii-1)} 'TE' polarization{indices(ii)+ii:end}}; % Proably error
        else
            polarization = {polarization{1:indices(ii)+(ii-1)} 'TM' polarization{indices(ii)+ii:end}};
        end
        nA = [nA(1:indices(ii)+(ii-1)) NaN nA(indices(ii)+ii:end)];
        kA = [kA(1:indices(ii)+(ii-1)) NaN kA(indices(ii)+ii:end)];
    end

    wl = unique(wl);
    sTM = strcmp(polarization, 'TM');
    sTE = strcmp(polarization, 'TE');
    TM = nA(sTM);
    aTM = kA(sTM);
    TE = nA(sTE);
    aTE = kA(sTE);
    SimRes = SimResults(wl, TM, aTM, hWG, hOrganic);
end

