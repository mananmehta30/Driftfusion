%% Initialise Driftfusio0n
par_memristor = pc('Input_files/memristor');


%% Find equilibrium solutions for different surface recombination rates
soleq_memristor = equilibrate(par_memristor);

 tpoints=(2*(Vmax-Vmin)/10*k_scan)+1;
 k_scan=0.1;
 Vmax=5;
 Vmin=-5;
 
 
 % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan, 10, tpoints);


