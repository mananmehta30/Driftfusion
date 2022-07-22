initialise_df
%% Load parameters and customize if necessary
% par = pc('Input_files/1_layer_test.csv');
par = pc('1_layer_MAPI_ITO_Ag.csv');

sc_array = [1e-12, 1e-4, 1e+4, 1e+12];
k_scan_array=[0.001;0.1];


for i = 1:length(sc_array)
    
    par.sc_r = sc_array(i);
    par.sc_l = par.sc_r;
    

    for j = 1:length(k_scan_array) %loop to run for different electrode workfunction
        
        
        par = refresh_device(par);      % This line is required to rebuild various arrays used DF
        
        %% Find equilibrium
        soleq(i, j) = equilibrate(par);
        dfplot.acx(soleq(i, j).ion)
        
        %% Current-voltage scan
        k_scan = k_scan_array(j);
        Vmax = 1.2;
        Vmin = -1.2;
        
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        sol_CV(i, j) = doCV(soleq(i, j).ion, 0, 0, Vmax, Vmin, k_scan, 1, 241);
 
    end   
end
