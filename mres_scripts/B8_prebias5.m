% Freeze ions at prebias updated version

% load parameters
par = pc('1_layer_MAPI_Ag_Al.csv');

% set high sc
sc = 1;

% sc array
sc_l_arr = sc*[0, 0, 1, 1];
sc_r_arr = sc*[0, 1, 0, 1];
% leg_txt = 

% set number of mobile ions
ai = 1;
par.N_ionic_species = ai;

% create logical string for anions on/off
boo_ani = string(par.N_ionic_species == 2);

% refresh device
par = refresh_device(par);

% Steady state CV input parameters
V0   =  0;
Vmax =  1;
Vmin =  -1;
cycles =1;
ppV = 100;

% set slow scan rate
k = 1e-8;

% adapt no. of points for equal spacing
points = (Vmax - Vmin)*2*ppV*cycles + 1;

% legend text
leg_txt = ["el only","no sc","sc_r = 1","sc_l = 1","sc_r = sc_l = 1"];

%% Approximate steady-state CV
% el only steady state
sol_eq = equilibrate(par);
sol_CV_SS(1) = doCV(sol_eq.el, 0, V0, Vmax, Vmin, k, cycles, points);

for sci = 2:5
   % Set BC
   par.sc_l = sc_l_arr(sci-1);
   par.sc_r = sc_r_arr(sci-1);
   
   par = refresh_device(par);
   
   sol_eq = equilibrate(par);
      
   
    % Aproximate steady-state CV
    sol_CV_SS(sci) = doCV(sol_eq.ion, 0, V0, Vmax, Vmin, k, cycles, points);

    % create diagram title
    title_SS = sprintf(['%i cycles, sc_l = %g, sc_r = %g, rate = %g,',...
        ' mobile anions: %s'], cycles, sol_CV_SS(sci).par.sc_l, ...
        sol_CV_SS(sci).par.sc_r, k, boo_ani);
end

%% plot steady-state CV
% plot CV
colors = prism(5).*.75;

figure
for sci = 1:length(sol_CV_SS)
    hold on
    dfplot.JtotVapp(sol_CV_SS(sci), 0)
end

plot([0.2 0.2],[1e-5 150],'--','color','#CCCCCC')

hold off
xlim([Vmin,Vmax])
legend(leg_txt)
title('steady state CV')
set(findall(gcf,'-property','FontSize'),'FontSize',25)
set(gca,'ColorOrder',colors,'yscale','log')
xlim([0 1])
plotbrowser

%% plot total ion concentrations
figure
for sci = 1:length(sol_CV_SS)
    hold on
    ACV(sol_CV_SS(sci));

    title('ion concentrations')

    plotbrowser
end
legend(repelem(leg_txt,2))


%% Freezing ions at set voltages

% freeze voltage
Vpre = 0.8;

for Vi = 1:length(sol_CV_SS)
    % Instead of recalculating every voltage, cut sol_CV_SS at desired points
    Vapp_SS = dfana.calcVapp(sol_CV_SS(Vi));

    % pick freeze voltage
    fprintf(['\nVpre = ',num2str(Vpre),' V\n'])

    % Get time points, forwards & backwards bias
    % nearest value
    [~, minind] = min(abs(Vapp_SS - Vpre));

    % corresponding time point
    tpre(Vi) = sol_CV_SS(Vi).t(minind);

    % extract single time point as prebias solution and freeze ions
    sol_pre(Vi) = extract_IC(sol_CV_SS(Vi), tpre(Vi));
    sol_pre(Vi).par.mobseti = 0;
    
    sol_CV_fr(Vi) = doCV(sol_pre(Vi), 0, V0, 1.2, -1.2, k, 1, 480);
end

sol_pre(Vi).par.mobseti = 1;

sound([0.2*sin(0.8*(1:500)), zeros(1,500), 0.2*sin(0.8*(1:500))],4000)


sol_CV_fr_pos = sol_CV_fr;
%% Plot all CVs

% exclude empty solution
succ = true(1, length(sol_CV_fr));

figure
for Vi = 1:length(sol_CV_fr)
    % plot CV
    hold on
    try
        dfplot.JtotVapp(sol_CV_fr(Vi), 0)
    catch
        warning('\nNo plot output at %g V \n', Vi);
        succ(Vi) = false;
    end
    
end

% Vsucc = Vpre_arr(succ ~= 0);
% leg_txt = string(Vsucc) + " V";

hold off
legend(leg_txt)
title(sprintf('Frozen ions at %g, mobile anions %s',[Vpre, boo_ani]));
plotbrowser

% sound([0.2*sin(0.8*(1:500)), zeros(1,500), 0.2*sin(0.8*(1:500))],4000)

%% stack figure overview
colors = repelem(prism(5),2,1).*.75;


figure
for iii = 1:5
% figure
% subplot(3,1,1)
% dfplot.ELx(sol_CV_fr(iii),0)
% set(gca,'Xticklabels',[]);
% a1 = gca;
% xlabel([])

subplot(2,1,1)
dfplot.npx(sol_CV_fr(iii),0)
set(gca,'Xticklabels',[], 'Yscale','log');
a2 = gca;
xlabel([])
lines1 = findobj(a2,'Type','line');
ylim([1e-9 1e21])
hold on

% a2.LineStyleorder{


subplot(2,1,2)
dfplot.acx(sol_CV_fr(iii),0)
a3 = gca;
set(gca, 'yscale','log')
lines2 = get(gca, 'Children');
ylim auto
hold on

linkaxes([a2,a3],'x')

% xlim([0 1])
% set([a1,a2,a3],'xscale','log')

% sgtitle(sprintf('V_{pre} = %g, %s',[Vpre,leg_txt(iii)]))

ha = get(gcf,'children');


end


hold off

delete(lines2([2,4,6,8,10]))
legend(['',leg_txt],'Location','eastoutside')
a3.ColorOrder = colors(1:2:9,:);


subplot(2,1,1)
legend(["",["n","p","n","p","n","p","n","p","n","p"] + ", " + repelem(leg_txt,2)],'Location','eastoutside')
a2.ColorOrder = colors;
set(lines1(1:2:9), 'LineStyle' ,'--');

% set(ha(7),'position',[.2 .7 .6 .3])
set(ha(4),'position',[.2 .53 .6 .45])
set(ha(2),'position',[.2 .08 .6 .45])

%% bulk cation concentration

figure
for ii = 1:5
dfplot.acx(sol_CV_fr(ii))
hold on
end

lines2 = get(gca, 'Children');
delete(lines2([2,4,6,8,10]))

legend(leg_txt)
xlim([100 300])
ylim auto

%% calculate conductance
sol_CV_fr = sol_CV_fr_pos(5);

[~,~,x,~,~,n] = dfana.splitsol(sol_CV_fr);

n = n(1,:);

% average
n_ave = sum(n)/length(n);

sig_ave = sol_CV_fr.par.e*sol_CV_fr.par.mue*n_ave;
GoA_ave = sig_ave/sol_CV_fr.par.d

% integral
sig_x = sol_CV_fr.par.e*sol_CV_fr.par.mue*n;
rho_x = 1./sig_x;

RoA_int = trapz(x, rho_x);
GoA_int = 1/RoA_int

GoA_int = (trapz(x, 1./sig_x))^-1


%% ELx at frozen pre-bias 0V
sol1 = sol_CV_fr(1);
sol2 = sol_CV_fr(2);
sol3 = sol_CV_fr(3);
sol4 = sol_CV_fr(4);
sol5 = sol_CV_fr(5);
t = 0.5/k;

figure
tl = tiledlayout(1,5);

nexttile
dfplot.ELx(sol1, t)
pbaspect([1 1 1])
% set(gca, 'xticklabels',[])
xlabel([])
a11 = gca;
title('el only')

nexttile
dfplot.ELx(sol2, t)
pbaspect([1 1 1])
% set(gca, 'xticklabels',[])
xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a12 = gca;
title('no sc')

nexttile
dfplot.ELx(sol3, t)
pbaspect([1 1 1])
% set(gca, 'xticklabels',[])
% xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a13 = gca;
title('sc_r')

nexttile
dfplot.ELx(sol4, t)
pbaspect([1 1 1])
legend('Location','eastoutside')
% set(gca, 'xticklabels',[])
xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a14 = gca;
s=findobj('type','legend');
delete(s)
title('sc_l')


nexttile
dfplot.ELx(sol5, t)
pbaspect([1 1 1])
legend('Location','eastoutside')
% set(gca, 'xticklabels',[])
xlabel([])
ylabel([])
set(gca, 'yticklabels',[])
a15 = gca;
title('sc_r , sc_l')


nexttile(5)
xx = legend('location','eastoutside');

tl.TileSpacing = 'compact';
tl.Padding = 'compact';
set(findall(gcf,'-property','FontSize'),'FontSize',30)

linkaxes([a11, a12, a13, a14, a15],'y')
% linkaxes([a21, a22, a23, a24, a25],'y')