function [Ex] = GenerateFieldSlot(model, order)
%G Summary of this function goes here
%   Detailed explanation goes here
    Params_NonSwept = mphgetexpressions(model.param); 
    Parameters = struct();

    for idx_param_list = 1:size(Params_NonSwept, 1)
        Parameters.(Params_NonSwept{idx_param_list, 1}) = struct('value', abs(Params_NonSwept{idx_param_list, 4}), 'phase', angle(Params_NonSwept{idx_param_list, 4}), 'unit', ['[' Params_NonSwept{idx_param_list, 5} ']']);
    end

    wSim = Parameters.wSim.value; 
    wWG = Parameters.wWG.value; 
    hSubstrate = Parameters.hSubstrate.value; 
    hRidge = Parameters.hRidge.value;
    hWG = Parameters.hWG.value; 
    hOrganic = Parameters.hOrganic.value; 
    hAir = Parameters.hAir.value; 
    dSlot = Parameters.dSlot.value; 
    n = 400; 
    x_edge = linspace(0, wSim, n); 
    y_edge = linspace(0, hSubstrate + hRidge + hWG + hOrganic + hAir, n); 
    [xx, yy] = meshgrid(x_edge, y_edge); 
    z = zeros(size(xx)); 
    for ii = 1:n
        for kk = 1:n
            ycond1 = yy(ii, kk) > hSubstrate & yy(ii, kk) < hSubstrate + hWG; 
            if ycond1
                xcond1 = xx(ii, kk) >= (wSim - dSlot)/2 - wWG & xx(ii, kk) < (wSim - dSlot)/2 - wWG + wWG/2;
                xcond2 = xx(ii, kk) >= (wSim - dSlot)/2 - wWG + wWG/2 & xx(ii, kk) < (wSim - dSlot)/2; 
                xcond3 = xx(ii, kk) >= (wSim - dSlot)/2 & xx(ii, kk) < (wSim + dSlot)/2; 
                xcond4 = xx(ii, kk) >= (wSim + dSlot)/2 & xx(ii, kk) < (wSim + dSlot)/2 + wWG/2; 
                xcond5 = xx(ii, kk) >= (wSim + dSlot)/2 + wWG/2 & xx(ii, kk) < (wSim + dSlot)/2 + wWG; 
                
                if xcond1
                    z(ii, kk) = -1; 
                elseif xcond2
                    z(ii, kk) = 1; 
                elseif xcond3
                    z(ii, kk) = 1; 
                elseif xcond4
                    z(ii, kk) = 1; 
                elseif xcond5
                    z(ii, kk) = -1; 
                end
            end
        end
    end
    Ex = z; 
end

