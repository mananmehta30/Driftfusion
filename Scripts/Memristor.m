initialise_df

%% Par file
par_memristor = pc('Input_files/memristor');
%% Equailibirum Solution
soleq_memristor = equilibrate(par_memristor);
%% Scan parameters
 Vmax=1.2;
 Vmin=-1.2;
k_scan=0.1;
%% Sol_CV
 % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)

% sol_CV = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan, 1, 100);
 sol_CV = doCV(soleq_memristor.el, 0, 0, Vmax, Vmin, k_scan, 1, 241);
%%
dfplot.JtotVapp(sol_CV,1);