% Load device and edit sc_r value
par = pc('1_layer_MAPI_Ag_Al.csv');

% scan rates for loops
% k_arr = sort([logspace(-3,0,4),5.*logspace(-3,0,4)]);
% k_arr = 1e-3;
% k_arr = [.01,.02,.03,.04,.05];
k_arr = 1e-1;

% sc_arr = logspace(-10,-6,5);
sc_arr = 1e-8;

% recalculate CV solutions
recalc = 1;

%% No sc case
%{


par_nosc = par;
par_nosc.sc_r = 0;
par_nosc = refresh_device(par_nosc);

sol_eq_nosc = equilibrate(par_nosc);


for ki = 1:length(k_arr)
    V0 = 0;
    Vmax = 1.2;
    Vmin = -1.2;
    cycles = 5;
    k = k_arr(ki);
    points = 480*cycles+1;
    
    sol_CV_nosc(ki) = doCV(sol_eq_nosc.ion, 0, V0, Vmax, Vmin, k, cycles, points);
    sol_CVel(ki) = doCV(sol_eq_nosc.el, 0, V0, Vmax, Vmin, k, 1, 481);
    
    
    figure
    dfplot.JtotVappC(sol_CV_nosc(ki),0)
    hold on 
    dfplot.JtotVapp(sol_CVel(ki),0)
    hold off
    title(sprintf('no sc, scan rate = %g Vs-1',k))
%     legend('No sc','el only')


    plotbrowser
end

%}

%% loop sc
for si = 1:length(sc_arr)
%     par.sc_l = sc_arr(si);
    par.sc_l = sc_arr(si);
    
    par = refresh_device(par);

    if recalc
        sol_eq = equilibrate(par);
    end
    
for ki = 1:length(k_arr)
% CV parameters
V0 = 0;
Vmax = 1.2;
Vmin = -1.2;
cycles = 3;
k = k_arr(ki);
points = 480*cycles+1;
% points = 5013;

if recalc
    % run el only 
    sol_CVel(si,ki) = doCV(sol_eq.el, 0, V0, Vmax, Vmin, k, 1, 481);

    % run full CV
    sol_CV(si,ki) = doCV(sol_eq.ion, 0, V0, Vmax, Vmin, k, cycles, points);

end

% split solution to get time array
[~,tarr] = dfana.splitsol(sol_CV(si,ki));

% % test time array against expected values (problem now fixed)
%{
tmax_target = ((Vmax-Vmin)*2*cycles)/k;
tmax_actual = max(tarr);

if abs(tmax_target-tmax_actual) > eps(tmax_target)
    warning(['Simulated max time does not match desired max time \n',...
        '-------- Simulated scan rate does not match desired scan rate \n',...
        '-------- \n','-------- \n','-------- \n','-------- \n'],'')
end
%}
% %{ 
%% plot JV curve
figure()

dfplot.JtotVappC(sol_CV(si,ki),0)
% hold on
% dfplot.JtotVapp(sol_CV,0)
hold on
dfplot.JtotVapp(sol_CVel(si,ki), 0)
hold off

title(sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1'...
    ,[cycles, par.sc_r, par.sc_l, k]))

plotbrowser

%% Find maximum of JV curve
J = dfana.calcJ(sol_CV(si,ki)).tot;
J = J(:,1);

[Jmax,Jind] = max(J);

Vapp = dfana.calcVapp(sol_CV(si,ki));
V_Jmax = Vapp(Jind);
% V_Jmax = 0.3;

t = V_Jmax/k;

%% carrier densities
% %{
figure()
subplot(1,2,1)
dfplot.acx(sol_CV(si,ki),t)
% title(sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1, Jmax at %f V',...
%     [cycles, par.sc_r, par.sc_l, k0, V_Jmax]))
legend('','anion','cation')
set(gca,'yscale','log')
%}

%% ELx
% %{
subplot(1,2,2)
dfplot.ELx(sol_CV(si,ki),t)
sgtitle(sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1, Jmax at %f V',...
    [cycles, par.sc_r, par.sc_l, k, V_Jmax]))

plotbrowser
%}

%% integrate ion concentrations over device
%%{
% input parameters
sol = sol_CV(si,ki);
cycles = cycles;

% plots total ion concentrations over time

% split solution to find necessary arrays
[~,t,x,~,~,n,p,a,c,~] = dfana.splitsol(sol);

Vapp = dfana.calcVapp(sol);
V0 = Vapp(1);
Vmax = max(Vapp);
Vmin = min(Vapp);

% create voltage axis tick labels
xtick_label_V = [repmat([V0, Vmax, V0, Vmin],1,cycles),V0];


for i = 1:length(c(:,1))
    integc(i) = trapz(x,c(i,:));
    intega(i) = trapz(x,a(i,:));
end

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

ax2 = axes(tl);
ax2.XAxisLocation = 'top';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
set(gca, 'YTick', []);
set(gca, 'xtick', [0:length(xtick_label_V)-1]./length(xtick_label_V));
set(gca, 'xticklabel', xtick_label_V)
xlabel('Voltage [V]')

title(sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1, V_{max} = %.2f V, V_{min} = %.2f V',...
    [cycles, par.sc_r, par.sc_l, k, Vmax, Vmin]))


plotbrowser
%}

end

end

