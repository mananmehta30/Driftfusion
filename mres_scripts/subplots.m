% needs B4_cycle_labelling to work
figure
t0 = tiledlayout(1,3);

nexttile
dfplot.JtotVappC(sol_CV(si,ki),0)
% hold on
% dfplot.JtotVapp(sol_CV,0)
hold on
dfplot.JtotVapp(sol_CVel(si,ki), 0)
hold off

% title(sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1'...
%     ,[cycles, par.sc_r, par.sc_l, k]))

xlim([-.3 1.2])

nexttile([1 2])

plot(t,integc)
hold on
plot(t,intega)
hold off

xlabel('time [s]')
ylabel('total number of ions [cm-2]')
legend('cations','anions')

set(gca, 'xtick', [0:12:144]);
set(gca, 'xticklabel', xtick_label_V)
xlabel('Voltage [V]')

% title(sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1, V_{max} = %.2f V, V_{min} = %.2f V',...
%     [cycles, par.sc_r, par.sc_l, k, Vmax, Vmin]))

xlim([0 144])

plotbrowser

%% B4 cycle labelling acx plots
figure
dfplot.acx(sol_CV, 25.7)
title('ions at -0.17 V, -6.1e-4 Acm^{-2}')
pbaspect([1 1 1])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
legend('location','eastoutside')
set(gca,'yscale','log')

figure
dfplot.acx(sol_CV, 46.4)
title('ions at -0.16 V, -5.9e-6 Acm^{-2}')
pbaspect([1 1 1])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
legend('location','eastoutside')
set(gca,'yscale','log')

%% B8 ELx and acx at 0.5 V steady state
sol1 = sol_CV_SS(1);
sol2 = sol_CV_SS(2);
sol3 = sol_CV_SS(3);
sol4 = sol_CV_SS(4);
sol5 = sol_CV_SS(5);
t = 0.5/k;

figure
tl = tiledlayout(2,5);

nexttile
dfplot.ELx(sol1, t)
pbaspect([1 1 1])
set(gca, 'xticklabels',[])
xlabel([])
a11 = gca;

nexttile
dfplot.ELx(sol2, t)
pbaspect([1 1 1])
set(gca, 'xticklabels',[])
xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a12 = gca;

nexttile
dfplot.ELx(sol3, t)
pbaspect([1 1 1])
legend('Location','eastoutside')
set(gca, 'xticklabels',[])
xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a13 = gca;

nexttile
dfplot.ELx(sol4, t)
pbaspect([1 1 1])
legend('Location','eastoutside')
set(gca, 'xticklabels',[])
xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a14 = gca;

nexttile
dfplot.ELx(sol5, t)
pbaspect([1 1 1])
legend('Location','eastoutside')
set(gca, 'xticklabels',[])
xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a15 = gca;

nexttile
dfplot.acx(sol1, t)

set(gca,'yscale','log')
ylim auto
pbaspect([1 1 1])
a21 = gca;

nexttile
dfplot.acx(sol2, t)
set(gca,'yscale','log')
ylim auto
pbaspect([1 1 1])
ylabel([])
set(gca, 'yticklabels',[])
a22 = gca;

nexttile
dfplot.acx(sol3, t)
set(gca,'yscale','log')
ylim auto
pbaspect([1 1 1])
ylabel([])
set(gca, 'yticklabels',[])
a23 = gca;

nexttile
dfplot.acx(sol4, t)
set(gca,'yscale','log')
ylim auto
pbaspect([1 1 1])
s=findobj('type','legend');
delete(s)
ylabel([])
set(gca, 'yticklabels',[])
a24 = gca;

nexttile
dfplot.acx(sol5, t)
set(gca,'yscale','log')
ylim auto
legend('','a','c','location','eastoutside')
pbaspect([1 1 1])
ylabel([])
set(gca, 'yticklabels',[])
a25 = gca;

nexttile(5)
xx = legend('location','eastoutside');

tl.TileSpacing = 'none';
tl.Padding = 'compact';
set(findall(gcf,'-property','FontSize'),'FontSize',30)

linkaxes([a11, a12, a13, a14, a15],'y')
linkaxes([a21, a22, a23, a24, a25],'y')