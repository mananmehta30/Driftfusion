% Simple memristortttt
initialise_df

%% Define memristor
par_memristor = pc('Input_files/memristor.csv');

%% Get Equilbrium solutions
soleq_memristor = equilibrate(par_memristor);


%% Cyclic Voltammogram scan
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
 k_scan = 1;
 tpoints=200;
 
Vmax = 5;
 Vmin = -5;
%%
sol_CV = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
%%
sol_CV_el = doCV(soleq_memristor.el, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);


%% Plot for ions
dfplot.JtotVapp(sol_CV, 0);
%set(gca,'YScale','log');

%% Plot for electrons
dfplot.JtotVapp(sol_CV_el, 0);
%set(gca,'YScale','log');