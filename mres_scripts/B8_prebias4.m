% Freeze ions at prebias updated version

% load parameters
par = pc('1_layer_MAPI_Ag_Al.csv');

% set high sc
sc = 1;
par.sc_l = sc;
par.sc_r = sc;

% set number of mobile ions


ai = 1;
par.N_ionic_species = ai;

% create logical string for anions on/off
boo_ani = string(par.N_ionic_species == 2);

% refresh device
par = refresh_device(par);

% initial equilibrium solution
sol_eq = equilibrate(par);


%% Aproximate steady-state CV

% input parameters
V0   =  0;
Vmax =  1;
Vmin =  -1;
cycles =1;
ppV = 100;

% set slow scan rate
k = 1e-8;

% adapt no. of points for equal spacing
points = (Vmax - Vmin)*2*ppV*cycles + 1;

% run CV
sol_CV_SS = doCV(sol_eq.ion, 0, V0, Vmax, Vmin, k, cycles, points);

% create diagram title
title_SS = sprintf(['%i cycles, sc_l = %g, sc_r = %g, rate = %g,',...
    ' mobile anions: %s'], cycles, sol_CV_SS.par.sc_l, ...
    sol_CV_SS.par.sc_r, k, boo_ani);

%% plot steady-state CV
% plot solution
figure
dfplot.JtotVapp(sol_CV_SS, 0)

title(title_SS)
xlim([Vmin,Vmax])

plotbrowser

% plot concentrations
ACt(sol_CV_SS, cycles, 0);

title(title_SS)

plotbrowser


%% Freezing ions at set voltages
% Instead of recalculating every voltage, cut sol_CV_SS at desired points

Vapp_SS = dfana.calcVapp(sol_CV_SS);

% array of all desired freeze voltage
Vpre_arr = linspace(Vmin, Vmax, abs(Vmax-Vmin)*5+1);
% Vpre_arr = 1;

% loop over all prebias voltages
for Vi = 1:length(Vpre_arr)

    % pick freeze voltage
    Vpre = Vpre_arr(Vi);
    fprintf(['\nVpre = ',num2str(Vpre),' V\n'])

    % Get time points, forwards & backwards bias
    % nearest value
    [~, minind] = min(abs(Vapp_SS - Vpre));

    % corresponding time point
    tpre(Vi) = sol_CV_SS.t(minind);

    % extract single time point as prebias solution and freeze ions
    sol_pre(Vi) = extract_IC(sol_CV_SS, tpre(Vi));
    sol_pre(Vi).par.mobseti = 0;
    
    sol_CV_fr(Vi) = doCV(sol_pre(Vi), 0, V0, 1.2, -1.2, k, 1, 480);

end

sol_pre(Vi).par.mobseti = 1;

sound([0.2*sin(0.8*(1:500)), zeros(1,500), 0.2*sin(0.8*(1:500))],4000)

%% Plot all CVs

% exclude empty solution
succ = true(1, length(Vpre_arr));

figure
for Vi = 1:length(Vpre_arr)
    % plot CV
    hold on
    try
        dfplot.JtotVapp(sol_CV_fr(Vi), 0)
    catch
        warning('\nNo plot output at %g V \n', Vpre_arr(Vi));
        succ(Vi) = false;
    end
    
end

Vsucc = Vpre_arr(succ ~= 0);
leg_txt = string(Vsucc) + " V";

hold off
legend(leg_txt)
title(sprintf('Frozen ions at voltages, mobile anions %s',boo_ani));
plotbrowser

% sound([0.2*sin(0.8*(1:500)), zeros(1,500), 0.2*sin(0.8*(1:500))],4000)

%% ELx, acx, npx of frozen ions
% diagrams with frozen ions at 0V after prebias

for di = 1: length(sol_CV_fr)
    
    figure
    subplot(1,3,1)
    dfplot.ELx(sol_CV_fr(1,di,:,:), 0)
    set(gca, 'xscale', 'log')
    
    subplot(1,3,2)
    dfplot.acx(sol_CV_fr(1,di,:,:), 0)
    set(gca, 'xscale', 'log')
    
    subplot(1,3,3)
    dfplot.npx(sol_CV_fr(1,di,:,:), 0)
    set(gca, 'xscale', 'log')    
    
    sgtitle(sprintf('freeze voltage %g V, mobile anions %s', Vpre_arr(di),...
        boo_ani))
   
    plotbrowser
end

%% conductance around 0V

for si = find(succ)
    
    Vapp_fr{si} = dfana.calcVapp(sol_CV_fr(si));
    J_arr = dfana.calcJ(sol_CV_fr(si));
    Jtot_fr{si} = J_arr.tot(:,1);
    % Get index of zero crossing
    indices = find(~any(Vapp_fr{si},1));
    ind(si) = indices(1);
    % Get JV gradient
    JV_grad = gradient(Jtot_fr{si},Vapp_fr{si});
    JV_grad0(si) = JV_grad(ind(si));


end

JV_grad0 = JV_grad0(succ);

figure
plot(Vsucc, abs(JV_grad0), 'x')
set(gca,'yscale', 'log')
title('Conductance at 0 V after applied pre-bias')
xlabel('Pre-bias [V]')
ylabel('Conductance [\Omega^{-1}cm^{-2}]')

%{
%% Conductance in Ohmic region
% loop over all freeze voltages
figure
for si = find(succ)
    Vpre = Vpre_arr(si);
    
    % onset voltage of linear region
    if Vpre >= 0
        Von = Vpre + 0.1;
    else
        Von = 0.5;
    end

    try
        Vfit{si} = Vapp_fr{si}(Vapp_fr{si} >= Von);
        Vfit{si} = Vfit{si}(1:floor(end/2)+1);
        Jfit{si} = Jtot_fr{si}(Vapp_fr{si} >= Von);
        Jfit{si} = Jfit{si}(1:floor(end/2)+1);
        
        l_fit{si} = polyfit(Vfit{si}, Jfit{si}, 1);

    catch
        warning(sprintf('no output for %g V', Vpre))
    end
    
    
    cond_ohm(si) = l_fit{si}(1);
    xfit = linspace(Von, 1.2);
    yfit = polyval(l_fit{si}, xfit);

    plot(xfit, yfit);
    hold on

end
hold off


sound([0.2*sin(0.8*(1:500)), zeros(1,500), 0.2*sin(0.8*(1:500))])
%}