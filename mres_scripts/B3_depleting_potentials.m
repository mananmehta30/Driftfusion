
% Load parameters
par = pc('1_layer_MAPI_ITO_Ag.csv');

% Calculate equilibrium solution
sol_eq = equilibrate(par);

%%
% Create list of right work function values to vary BI potential
% Phi_right_arr = [-4.0, -4.25, -4.5, -4.75, -5.0];
Phi_right_arr = -4.8;

% Temporary parameters variable
par_temp = par;
% par_temp.N_ionic_species = 2;

% Change WF and refresh
par_temp.Phi_right = Phi_right_arr;
par_temp = refresh_device(par_temp);

% Initial solution and CV
sol_eq = equilibrate(par_temp);
sol_CV = doCV(sol_eq.ion, 0, 0, 1.2, -1.2, 1e-6, 2, 481);

% Plot
figure()
sgtitle(['\Delta\phi = ',num2str(abs(par_temp.Phi_left-par_temp.Phi_right)),...
    ' V, No. mobile ions: ',num2str(par_temp.N_ionic_species)])

% JV Curve
subplot(1,3,1)
dfplot.JtotVapp(sol_CV,0)

% ELx at 0 V
subplot(1,3,2)
dfplot.ELx(sol_CV, 0)

% ion concentrations, 0 V
subplot(1,3,3)
dfplot.acx(sol_CV, 0)

