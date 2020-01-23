% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [] = PlotDispersionRelation(results, varargin)
%PLOTDISPERSIONRELATION Summary of this function goes here
%   Detailed explanation goes here
    % Collect and Convert Data
    options = struct(...
        'fileName', '');
    optionNames=fieldnames(options);

    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
        if any(strcmp(inpName,optionNames))
            options(1).(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
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


    % Plot da shit
    % Dispersion Relation
    figure(1)
    hold on;
    plot(1.24./(wl*1e6), TE, 'LineWidth', 2, 'Color', [1, 0.643, 0.643]);
    plot(1.24./(wl*1e6), TM, 'LineWidth', 2, 'Color', [0.655, 0, 0]);
    hold off;
    xlabel('Photon Energy [eV]', 'FontSize', 14)
    ylabel('Effective Refractive Index', 'FontSize', 14)
    legend('TE00', 'TM00')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)
    ax = gca;
    labels = string(ax.XAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.XAxis.TickLabels = labels;
    labels = string(ax.YAxis.TickLabels);
    labels(1:2:end) = nan;
    ax.YAxis.TickLabels = labels;
    if save
        saveas(gcf, ['./Figures/' options.('fileName') '_RefractiveIndex.svg'])
        saveas(gcf, ['./Figures/' options.('fileName') '_RefractiveIndex.png'])
    end
    % Calculate Ideal Curve
    k0 = TE(1)*2*pi/wl(1)-TE(1)*2*pi/(2*wl(1));
    y = TE(1)*pi./wl(1)*ones(1, length(wl));
    % Momentum
    figure(2)
    hold on;
    plot(1.24./(wl*1e6), 1e-6*TE*2*pi./wl, 'LineWidth', 2, 'Color', [1, 0.643, 0.643]);
    plot(1.24./(wl*1e6), 1e-6*TM*2*pi./wl, 'LineWidth', 2, 'Color', [0.655, 0, 0]);
    %plot(1.24./(wl*1e6), 1e-6*y, 'LineWidth', 2, 'Color', [0.7969, 0, 0], 'LineStyle', '-');
    hold off;
    xlabel('Photon Energy [eV]', 'FontSize', 14)
    ylabel('Momentum [1/um]', 'FontSize', 14)
    legend('TE00', 'TM00')
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)
    ax = gca;
    labels = string(ax.XAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.XAxis.TickLabels = labels;
    labels = string(ax.YAxis.TickLabels);
    labels(1:2:end) = nan;
    ax.YAxis.TickLabels = labels;
    if save
        saveas(gcf, ['./Figures/' options.('fileName') '_Momentum.svg'])
        saveas(gcf, ['./Figures/' options.('fileName') '_Momentum.png'])
    end

    % Losses
    figure(3)
    hold on;
    plot(1.24./(wl*1e6), aTE*1e-2, 'LineWidth', 2, 'Color', [1, 0.643, 0.643]);
    plot(1.24./(wl*1e6), aTM*1e-2, 'LineWidth', 2, 'Color', [0.655, 0, 0]);
    hold off;
    xlabel('Photon Energy [eV]', 'FontSize', 12)
    ylabel('Losses [dB/cm]', 'FontSize', 12)
    legend('TE00', 'TM00')

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)

    ax = gca;
    labels = string(ax.XAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.XAxis.TickLabels = labels;
    labels = string(ax.YAxis.TickLabels);
    labels(1:2:end) = nan;
    ax.YAxis.TickLabels = labels;
    if save
        saveas(gcf, ['./Figures/' options.('fileName') '_Losses.svg'])
        saveas(gcf, ['./Figures/' options.('fileName') '_Losses.png'])
    end

    % Surface Plot for Momentum mismatch TM
    [ww, vv] = meshgrid(wl, wl);
    kp = TM*2*pi./ww;
    ks = TM*2*pi./vv;
    ki = TM*2*pi./vv;

    figure(4)
    contourf(1.24./(ww*1e6), 1.24./(vv*1e6), (kp-(ks + ki))*1e-6)
    c = colorbar();
    xlabel('Photon Energy (pump) [eV]')
    ylabel('Photon Energy (i, s) [eV]')

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)

    ax = gca;
    labels = string(ax.XAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.XAxis.TickLabels = labels;
    labels = string(ax.YAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.YAxis.TickLabels = labels;
    ylabel(c, 'Momentum Mismatch [1/um]')
    c.FontSize = 12;
    labels = string(c.TickLabels);
    labels(1:2:end) = nan;
    c.TickLabels = labels;
    if save
        print(['./Figures/' options.('fileName') '_MomentumMismatch'], '-dtiffn', '-r300')
        saveas(gcf, ['./Figures/' options.('fileName') '_MomentumMismatch.png'])
    end

    % Calculate Dispersion parameters
    omega = 2*pi*c_const./wl;
    dw = omega(1) - omega(2);
    beta1TM = zeros(1, length(omega)-2);
    beta1TE = zeros(1, length(omega)-2);
    beta2TM = zeros(1, length(omega)-2);
    beta2TE = zeros(1, length(omega)-2);
    for ii = 2:length(omega)-1
       beta1TM(ii-1) = 1/c_const*(TM(ii) + omega(ii)*(TM(ii) - TM(ii -1))/dw);
       beta1TE(ii-1) = 1/c_const*(TE(ii) + omega(ii)*(TE(ii) - TE(ii -1))/dw);
       beta2TM(ii-1) = 1/c_const*(2*(TM(ii + 1) - TM(ii -1))/dw + omega(ii)*(TM(ii + 1) - 2*TM(ii) + TM(ii -1))/dw^2);
       beta2TE(ii-1) = 1/c_const*(2*(TE(ii + 1) - TE(ii -1))/dw + omega(ii)*(TE(ii + 1) - 2*TE(ii) + TE(ii -1))/dw^2);
    end

    figure(5)
    hold on;
    plot(omega(2:end-1), 1./beta1TE, 'LineWidth', 2, 'Color', [1, 0.643, 0.643]);
    plot(omega(2:end-1), 1./beta1TM, 'LineWidth', 2, 'Color', [0.655, 0, 0]);
    hold off;
    xlabel('Photon Energy [eV]', 'FontSize', 12)
    ylabel('Dispersion  [1e8 m/s]', 'FontSize', 12)
    legend('TE00', 'TM00')

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)

    ax = gca;
    labels = string(ax.XAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.XAxis.TickLabels = labels;
    labels = string(ax.YAxis.TickLabels);
    labels(1:2:end) = nan;
    ax.YAxis.TickLabels = labels;
    if save
        saveas(gcf, ['./Figures/' options.('fileName') '_GV.svg'])
        saveas(gcf, ['./Figures/' options.('fileName') '_GV.png'])
    end

    figure(6)
    hold on;
    plot(omega(2:end-1), -1e6*omega(2:end-1).^2/(2*pi*c_const).*beta2TE, 'LineWidth', 2, 'Color', [1, 0.643, 0.643]);
    plot(omega(2:end-1), -1e6*omega(2:end-1).^2/(2*pi*c_const).*beta2TM, 'LineWidth', 2, 'Color', [0.655, 0, 0]);
    hold off;
    xlabel('Photon Energy [eV]', 'FontSize', 12)
    ylabel('Group Velocity Dispersion [ps/(nm km)]', 'FontSize', 12)
    legend('TE00', 'TM00')

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)

    ax = gca;
    labels = string(ax.XAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.XAxis.TickLabels = labels;
    labels = string(ax.YAxis.TickLabels);
    labels(1:2:end) = nan;
    ax.YAxis.TickLabels = labels;
    if save
        saveas(gcf, ['./Figures/' options.('fileName') '_GVD.svg'])
        saveas(gcf, ['./Figures/' options.('fileName') '_GVD.png'])
    end

    figure(7)
    contourf(1.24./(ww*1e6), 1.24./(vv*1e6), 1e6*pi./abs((kp-(ks + ki))))
    c = colorbar();
    xlabel('Photon Energy (pump) [eV]')
    ylabel('Photon Energy (i, s) [eV]')

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)

    ax = gca;
    labels = string(ax.XAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.XAxis.TickLabels = labels;
    labels = string(ax.YAxis.TickLabels);
    labels(2:2:end) = nan;
    ax.YAxis.TickLabels = labels;

    c.FontSize = 12;
    set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%0.3f')) )
%     labels = string(c.TickLabels);
%     labels(1:2:end) = nan;
%     c.TickLabels = labels;
    ylabel(c, 'Coherence Length [um]');
    if save
        print(['./Figures/' options.('fileName') '_CoherenceLength'], '-dtiffn', '-r300')
        saveas(gcf, ['./Figures/' options.('fileName') '_MomentumMismatch.png'])
    end
    % Prepare Data for scatter3 plot

    [wlx, wly, wlz] = meshgrid(wl, wl, wl);
    mismatch = zeros(length(wl), length(wl), length(wl));
    for ii = 1:length(wl)
        for kk = 1:length(wl)
            for ll = 1:length(wl)
                dk = TM(ii)*2*pi/wl(ii) - (TM(kk)*2*pi/wl(kk) + TM(ll)*2*pi/wl(ll));
                if abs(dk) < 2e6
                    mismatch(ii, kk, ll) = dk;
                else
                    mismatch(ii, kk, ll) = NaN;
                end
            end
        end
    end
    figure(8)

    figure()
    scatter3(wlx(:)*1e6, wly(:)*1e6, wlz(:)*1e6, 30, mismatch(:)*1e-6, 'filled');
    hold on
    wln = min(wly(:))*ones(size(wly));
    wlm = max(wly(:))*ones(size(wly)) + 0.3*1e-6;
    scatter3(wlx(:)*1e6, wly(:)*1e6, wln(:)*1e6, 30, mismatch(:)*1e-6, 'filled');
    scatter3(wlx(:)*1e6, wlm(:)*1e6, wlz(:)*1e6, 30, mismatch(:)*1e-6, 'filled');
    scatter3(wlm(:)*1e6, wly(:)*1e7, wlz(:)*1e6, 30, mismatch(:)*1e-6, 'filled');
    xlim([min(wl)*1e6, max(wl)*1e6 + 0.300])
    ylim([min(wl)*1e6, max(wl)*1e6 + 0.300])
    zlim([min(wl)*1e6, max(wl)*1e6 + 0.300])
    hold off
    cbar = colorbar;
    ylabel('\lambda_s [um]', 'fontsize', 12)
    set(get(gca,'ylabel'),'rotation',-30)
    xlabel('\lambda_p [um]', 'fontsize', 12)
    set(get(gca,'xlabel'),'rotation',15)
    zlabel('\lambda_i [um]', 'fontsize', 12)
    cbar.Label.String= 'Phase Mismatch [um^{-1}]';
    cbar.Label.FontSize = 12;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)

    if save
        print(['./Figures/' options.('fileName') '_MomentumMismatch3D'], '-dtiffn', '-r300')
        saveas(gcf, ['./Figures/' options.('fileName') '_MomentumMismatch3D.png'])
    end
end

