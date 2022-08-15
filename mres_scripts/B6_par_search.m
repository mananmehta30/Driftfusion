% import parameters
par_ex = pc('1_layer_MAPI_ITO_Ag.csv');

par_ex.Phi_left = -5;

par_ex = refresh_device(par_ex);

%% parameter search (doesn't work yet)
% scan rate
% sc_r
% 
% sol_ex_k_scr = explore.explore2par(par_ex, {'sc_r','sc_l'},...
%     {logspace(-20,-5,4),logspace(-20,-5,4)},200);

%% JV 
% equilibrate device
sol_eq = equilibrate(par_ex);

% parameters
k = 1e-1;
points = 201;
Vstart = 0; 
Vend = 1.2;


sol_JV = doJV(sol_eq.ion, k, points, 0, 1, Vstart, Vend, 1);

% %{ 
%% plot
figure
dfplot.JtotVapp(sol_JV.dk.f,0)
hold on
dfplot.JtotVapp(sol_JV.dk.r,0)
hold off
%}