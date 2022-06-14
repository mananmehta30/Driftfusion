% Single-layer MaPbI device, variable no. of ionic carriers etc.
% 14/04/2021
%
%% Load parameters
par = pc('1_layer_MAPI_ITO_Ag.csv');

sc = 1e-4;
k = 1e-6;

par.sc_l = sc;
par.sc_r = sc;

par = refresh_device(par);

%% Equilibrium solution
soleq = equilibrate(par);

sol_CVel = doCV(soleq.el,0,0,1.2,-1.2,k,1,241);
sol_CV = doCV(soleq.ion,0,0,1.2,-1.2,k,2,481);

figure
dfplot.JtotVappC(sol_CV,0)
hold on
dfplot.JtotVapp(sol_CVel,0)
hold off