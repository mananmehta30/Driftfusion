
%%Initialize
initialise_df
%%
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

%% Plot distribution
k_scan_index=2;
for i = 1:length(sc_array)
   dfplot.acx(sol_CV(i,k_scan_index));
     legstr_acx{i} = ['SC rate =', num2str(sc_array(i)),'velocity units'];
     hold on
end
hold off
legend(legstr_acx)
title('Scan rate =',num2str(k_scan_array(k_scan_index)),'V/s');
%% Plot distribution
sc_index=2;
for i = 1:length(k_scan_array)
   dfplot.acx(sol_CV(sc_index,i));
     legstr_ac2x{i} = ['SC rate =', num2str(sc_array(i)),'V/s'];
     hold on
end
hold off
legend(legstr_ac2x)
title('SCrate =',num2str(sc_array(sc_index)),'velocity units');