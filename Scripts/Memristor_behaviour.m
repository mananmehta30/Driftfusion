%% Simulate Ag and Iodine reaction by varying vacancy concentration at the interface and introduce
% a surface recombination velocity
initialise_df

%% Define memristor
par_memristor = pc('Input_files/memristor.csv');

%% Get Equilbrium solutions
soleq_memristor = equilibrate(par_memristor);


%% Cyclic Voltammogram scan
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
 k_scan = 0.1;
Vmax = 1;
 Vmin = -1;
tpoints=(2*(Vmax-Vmin)/(10*k_scan))+1;
sol_CV = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);

%%No change with sol_JV achieved. Need to find how to vary the surface
%%recombination formula

