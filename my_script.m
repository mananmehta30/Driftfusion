%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Windows
% par_sio2 = pc('.\Input_files\pog2');
% Filepath Mac
par_alox = pc('Input_files/alox.csv');
%% Equilibrium solutions
soleq_alox = equilibrate(par_alox);

%% Plot equilibrium energy level diagram
dfplot.ELnpx(soleq_alox.ion)

%% Current-voltage scan
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% JVsol = doJV(soleq_sio2.ion, 100e-3, 201, 1, 0, 0, 1, 1);

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
k_scan = 0.1;
sol_CV = doCV(soleq_alox.ion, 0, 0, 1, -1, k_scan, 2, 201);

%% Plot JV scan
dfplot.JtotVapp(sol_CV, 0);
%set(gca,'YScale','log')

%% Plot anion and cation densities
dfplot.acx(sol_CV, 1/k_scan*[0, 0.5, 1.0, 2.5, 3.0]);
% 
% ylim([-30e-3,10e-3])
% xlim([-0.2, 1.2])

%% Plot Vapp vs time
dfplot.Vappt(sol_CV)
