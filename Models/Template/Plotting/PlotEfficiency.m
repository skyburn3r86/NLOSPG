% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [] = PlotEfficiency(x, y, nu)
%PLOTDISPERSIONRELATION Summary of this function goes here
%   Detailed explanation goes here
    % Collect and Convert Data
    save = true;
    [n, m] = size(nu);
    N =  length(nu(1, 1));
    hWG = x;
    hOrganic = y;
    [xx, yy] = meshgrid(hWG, hOrganic);
    % Plot da shit
    % Dispersion Relation
    for ii = 1:N
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
            saveas(gcf, ['./Figures/' options.('fileName') '_Losses_' num2str(wl(ii)) '_' num2str(ii) 'nm.svg'])
            saveas(gcf, ['./Figures/' options.('fileName') '_Losses_' num2str(wl(ii)) '_' num2str(ii) 'nm.png'])
        end
    end

end