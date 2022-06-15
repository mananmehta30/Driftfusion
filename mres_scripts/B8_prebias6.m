% equilibrium ELx for all boundary cases
%% Initialize Driftfusion
initialise_df

%%
par = pc('1_layer_MAPI_Ag_Al.csv');

% set high sc
sc = 1;

% sc array
sc_l_arr = sc*[0, 0, 1, 1];
sc_r_arr = sc*[0, 1, 0, 1]; 

leg_txt = ["el only","no sc","sc_r = 1","sc_l = 1","sc_r = sc_l = 1"];

%el only case
soleq(1) = equilibrate(par);

% no sc case
par.sc_r = 0;
par.sc_l = 0;
soleq(2) = soleq(1);

% scl case
par.sc_r = 0;
par.sc_l = 1;
par = refresh_device(par);

soleq(3) = equilibrate(par);

% scr case
par.sc_r = 1;
par.sc_l = 0;
par = refresh_device(par);

soleq(4) = equilibrate(par);

% scr,scl case
par.sc_r = 1;
par.sc_l = 1;
par = refresh_device(par);

soleq(5) = equilibrate(par);

%% plot ELx

figure
tl = tiledlayout(2,3);

nexttile
dfplot.ELx(soleq(1).el, 0)
pbaspect([1 1 1])
set(gca, 'xticklabels',[])
xlabel([])
% set(get(gcf,'children'),'Position',[0 0 1 1 ])

nexttile
dfplot.ELx(soleq(4).ion, 0)
pbaspect([1 1 1])
set(gca, 'xticklabels',[])
xlabel([])

nexttile
dfplot.ELx(soleq(5).ion, 0)
pbaspect([1 1 1])
legend('Location','eastoutside')
set(gca, 'xticklabels',[])
xlabel([])

nexttile
dfplot.acx(soleq(1).el,0)
hold on
dfplot.npx(soleq(1).el,0)
hold off
ylim auto
pbaspect([1 1 1])

nexttile
dfplot.acx(soleq(4).ion,0)
hold on
dfplot.npx(soleq(4).ion,0)
hold off
ylim auto
pbaspect([1 1 1])
s=findobj('type','legend');
delete(s)

nexttile
dfplot.acx(soleq(5).ion,0)
hold on
dfplot.npx(soleq(5).ion,0)
hold off
ylim auto
legend('','a','c','n','p','location','eastoutside')
pbaspect([1 1 1])

nexttile(3)
xx = legend('location','eastoutside');

tl.TileSpacing = 'none';
tl.Padding = 'compact';
set(findall(gcf,'-property','FontSize'),'FontSize',30)