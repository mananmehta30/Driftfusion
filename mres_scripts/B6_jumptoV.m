% Load device and edit sc_r value
par = pc('1_layer_MAPI_ITO_Ag.csv');

% Open Ag boundary to cations
par.sc_r = 1;

par = refresh_device(par);
% initial solution
sol_eq = equilibrate(par);

%% standard JV

% parameters
k = 1e-10;
Vstart = 0;
Vend = 5;
points = abs(Vend-Vstart)*200+1;
% sol_JV0 = doJV(sol_eq.ion, k, points, 0, 1, Vstart, Vend, 1);


%% plot
%{
figure
dfplot.JV(sol_JV0, 1)
title(sprintf(['sc_r = %.0e cms-1, sc_l = %.0e cms-1, ',...
    'rate = %.0e Vs-1'], [par.sc_r, par.sc_l, k]))
%}

%% jump to voltage
% set prebias voltage
V_pre = -2;
t_dwell = 500;


% for i = 1:length(loopvar)
    
sol_pre = jumptoV(sol_eq.ion, V_pre, t_dwell, 1, 0, 0, 0);

[~,t,x,~,~,n,p,a,c,~] = dfana.splitsol(sol_pre);

for i = 1:length(c(:,1))
    integc(i) = trapz(x,c(i,:));
    intega(i) = trapz(x,a(i,:));
end

figure
plot(t,integc)
hold on
plot(t,intega)
hold off

xlabel('time [s]')
ylabel('total number of ions [cm-2]')
legend('cations','anions')

%% Move voltage back to 0V
sol_backto0 = doJV(sol_pre(end,:,:), 1, abs(V_pre) * 100, 0, 0, V_pre, 0, 1);

%% Run JV with frozen ions
sol_JV = doJV(sol_backto0.dk.f, 1, points, 0, 0, 0, Vend, 1);
%%
figure
dfplot.JV(sol_JV, 1)