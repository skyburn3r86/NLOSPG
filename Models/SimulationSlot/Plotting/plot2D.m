% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [] = plot2D(x, y, z, varargin)
% Creates a ContourPlot
%   Detailed Summary Here
    options = struct(...
            'filename', '', ...
            'xlabel', '', ...
            'ylabel', '', ...
            'title', '');
        optionNames=fieldnames(options);

        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('Arguments needs propertyName/propertyValue pairs')
        end

        for pair = reshape(varargin,2,[])
            inpName = pair{1};
           if any(strcmp(inpName,optionNames))
              options(1).(inpName) = pair{2}{1};
              options(2).(inpName) = pair{2}{2};
           else
              error('%s is not a recognized parameter name',inpName)
           end
        end

    close all
    [xx, yy] = meshgrid(x, y);

    figure(1)
    pcolor(xx,yy,transpose(z))
    colormap
    xlabel(options.('xlabel'), 'FontSize', 12)
    ylabel(options.('ylabel'), 'FontSize', 12)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)

    if ~strcmp(options.('filename'), '')
        saveas(gcf, options.('filename'))
    end
end