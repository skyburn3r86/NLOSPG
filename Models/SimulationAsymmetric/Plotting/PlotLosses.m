% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [] = PlotLosses(losses)
%PLOTDISPERSIONRELATION Summary of this function goes here
%   Detailed explanation goes here
    % Collect and Convert Data
    save = true;
    [n, m] = size(losses);
    wl =  losses{1,1}.('wl');
    hWG = zeros(n*m, length(wl));
    hOrganic = zeros(n*m, length(wl));
    alpha = zeros(n*m, length(wl));
    for ii = 1:n
        for kk = 1:m
            alpha((ii-1)*m + kk, :) = losses{ii, kk}.('alpha');
            hOrganic((ii-1)*m + kk, :) = losses{ii, kk}.('hOrganic');
            hWG((ii-1)*m + kk, :) = losses{ii, kk}.('hWG');
        end
    end

    % Plot da shit
    % Dispersion Relation
    for ii = 1:length(wl)
        [xx, yy] = meshgrid(unique(hOrganic(:, ii)), unique(hWG(:, ii)));
        Z = griddata(hOrganic(:, ii), hWG(:, ii), alpha(:, ii), xx, yy);
        figure(ii)
        pcolor(xx,yy,Z)
        colormap hot
        shading interp
        xlabel('hWG')
        ylabel('hOrganic')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12)
        a = get(gca,'YTickLabel');
        set(gca,'YTickLabel',a,'fontsize',12)
        set(gca,'ytick',linspace(0, max(hWG), 5), 'fontsize', 12);
        set(gca,'xtick',linspace(0, max(hOrganic), 5), 'fontsize', 12);
        if save
            saveas(gcf, ['./Figures/' options.('fileName') '_Losses_' num2str(wl(ii)) 'nm.svg'])
            saveas(gcf, ['./Figures/' options.('fileName') '_Losses_' num2str(wl(ii)) 'nm.png'])
        end
    end
end

