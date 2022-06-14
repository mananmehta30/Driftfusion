function [A,C] = ACt(sol, cycles, V_ax)
    % plots total ion concentrations over time
    % outputs concentration arrays over time

    % split solution to find necessary arrays
    [~,t,x,~,~,~,~,a,c,~] = dfana.splitsol(sol);
    
    % get voltage points
    Vapp = dfana.calcVapp(sol);
    V0 = Vapp(1);
    Vmax = max(Vapp);
    Vmin = min(Vapp);

    % create voltage axis tick labels
    xtick_label_V = [repmat([V0, Vmax, V0, Vmin],1,cycles),V0];

    integc = nan(1,length(c(:,1)));
    intega = nan(1,length(a(:,1)));
    for i = 1:length(c(:,1))
        integc(i) = trapz(x,c(i,:));
        intega(i) = trapz(x,a(i,:));
    end

    A = intega;
    C = integc;
    
    figure
    tl = tiledlayout(1,1);
    ax1 = axes(tl);

    plot(t,integc)
    hold on
    plot(t,intega)
    hold off

    xlabel('time [s]')
    ylabel('total number of ions [cm-2]')
    legend('cations','anions')
    set(gca,'yscale','log')

    if V_ax
        ax2 = axes(tl);
        ax2.XAxisLocation = 'top';
        ax2.Color = 'none';
        ax1.Box = 'off';
        ax2.Box = 'off';
        set(gca, 'YTick', []);
        set(gca, 'xtick', [0:(length(xtick_label_V)-1)]./(length(xtick_label_V)-1));
        set(gca, 'xticklabel', xtick_label_V)
        xlabel('Voltage [V]')
    end
end