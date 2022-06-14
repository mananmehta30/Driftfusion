function ACV(sol)
    % plots total ion concentrations over time
    % outputs concentration arrays over time

    % split solution to find necessary arrays
    [~,~,x,~,~,~,~,a,c,~] = dfana.splitsol(sol);
    
    Vapp = dfana.calcVapp(sol);

    integc = nan(1,length(c(:,1)));
    intega = nan(1,length(a(:,1)));
    for i = 1:length(c(:,1))
        integc(i) = trapz(x,c(i,:));
        intega(i) = trapz(x,a(i,:));
    end

    A = intega;
    C = integc;
    
    plot(Vapp,integc)
    hold on
    plot(Vapp,intega)
    hold off

    xlabel('Voltage [V]')
    ylabel('total number of ions [cm-2]')
    legend('cations','anions')
    set(gca,'yscale','log')
end