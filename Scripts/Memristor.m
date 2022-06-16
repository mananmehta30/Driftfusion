%% Initialise Driftfusion

initialise_df;

%% Get file parameters
par_memristor = pc('Input_files/memristor');


%% Input Parameters

k_scan_array=[0.01,0.1,1];
sc_array = [0,1e-12, 1e-10, 1e-8, 1e-6, 1e-4];

 Vmax=5;
 Vmin=-5;
 
tpoints=241; 

%% Find equilibrium solutions for different surface recombination rates

% sn is Electron surface recombination velocity
% For electrode layers, entries for sn and sp are stored in the parameters object (par)
% as the distinct properties sn_l, sn_r, sp_l, and sp_r rather than as part
% of the sn and sp arrays
soleq_memristor = equilibrate(par_memristor);  
  

  %%
for i = 1:length(sc_array) % Loop to run for different recombination velocities
    par_memristor = refresh_device(par_memristor);
     df.sc_r(1) = sc_array(i);
    soleq_memristor(i) = equilibrate(par_memristor);   
end


%% Calculate electron only solution
el_CV = doCV(soleq(1).el, 0, 0, 1.2, -1.2, 1e-1, 2, 241);


%% Solve for different BC, both sides, medium scan rate (0.1 Vs-1), two cycles

% Graph shows that higher surface recombintion velocity i.e. higher
% recombintion rate leads to more hystersis and less current
k_scan_array = 0.1;
cycles = 1;

figure()
for i = 1:length(sc_array)
     % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
    sol_CV(i) = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan_array, 10, tpoints);
    dfplot.JtotVapp(sol_CV(i),0)
    hold on
end
%Total Current Plot for only electron
dfplot.JtotVapp(el_CV,0)
hold off
% set(gca,'yscale','log')
legentries = cellstr(num2str(sc_array', 'sc=%g'));
legentries{end+1} = 'el only';
legend(legentries)
title(sprintf('%i cycle, scan rate = %g Vs-1, sc on both sides',[cycles,k_scan_array]))


