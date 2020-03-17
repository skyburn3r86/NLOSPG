% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [x, y] = plot1D(x, y)
% Summary Here
%   Detailed Summary Here
    options = struct(...
                 'filename', '', ...
                 'xlabel', '', ...
                 'ylabel', '', ...
                 'title', '', ...
                 'plot', 'normal');

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
    figure(1)
    if strcmp(options.('plot'), 'semilogy')
        semilogy(x, y, 'LineWidth', 2.0, 'Color', [0.8, 0.0, 0.0])
    elseif strcmp(options.('plot'), 'semilogx')
        semilogx(x, y, 'LineWidth', 2.0, 'Color', [0.8, 0.0, 0.0])
    elseif strcmp(options.('plot'), 'loglog')
        loglog(x, abs(y), 'LineWidth', 2.0, 'Color', [0.8, 0.0, 0.0])
    else
        plot(x, y, 'LineWidth', 2.0, 'Color', [0.8, 0.0, 0.0])
    end

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