
par_memristor = pc('Input_files/memristor');

soleq_memristor = equilibrate(par_memristor);

 Vmax=5;
 Vmin=-5;
k_scan=0.1;

%%
 % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)

 sol_CV = doCV(soleq_memristor.el, 0, 0, Vmax, Vmin, k_scan, 10, 400);
%%
dfplot.JtotVapp(sol_CV,1);