%% Simulate Ag and Iodine reaction by varying vacancy concentration at the interface and introduce
% a surface recombination velocity
initialise_df

%% Define memristor
par_memristor = pc('Input_files/memristor.csv');
par_memristor_interface_reactions = par_memristor;
%% Get Equilbrium solutions
soleq_memristor = equilibrate(par_memristor);
soleq_memristor_interface_reactions = equilibrate(par_memristor_interface_reactions);


%% Cyclic Voltamogram scan
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
 k_scan = 0.1;
Vmax = 1;
 Vmin = -1;
tpoints=(2*(Vmax-Vmin)/(10*k_scan))+1;
sol_CV = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
dfplot.ELnpx(sol_CV, 0);

