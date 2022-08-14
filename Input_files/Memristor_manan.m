% Simple memristor
initialise_df
%% Define memristor Ag/Au
%par_memristor = pc('Input_files/memristor_silver_both_sides_400nm.csv');%Ag both sides

%par_memristor = pc('Input_files/memristor_gold_both_sides_400nm.csv');% Au both sides

%par_memristor = pc('Input_files/memristor_silver_left_gold_right.csv');% Ag left Au right side

par_memristor = pc('Input_files/memristor_gold_left_silver_right.csv');% Au left Ag right side

par_memristor.N_ionic_species=1;
%% Get Equilbrium solutions
soleq_memristor = equilibrate(par_memristor);

%% Equlibrium Plots
%dfplot.acx(soleq_memristor.ion);
%% Cyclic Voltammogram scan
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
 k_scan = 0.0001;
 tpoints=200;

Vmax = 0.5;
Vmin = -0.5;
Vmaxtarray=Vmax/k_scan;
Vmintarray=3.*Vmaxtarray;
%% Sol_CV
sol_CV = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);


%%
%sol_CV_el = doCV(soleq_memristor.el, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);


%% Plot for ions
dfplot.JtotVapp(sol_CV, 0);
%%
%set(gca,'YScale','log');
% dfplot.JVapp(sol_CV, 0)
%% Get direction
%  J = dfana.calcJ(sol_CV, "sub");
%  Vapp = dfana.calcVapp(sol_CV);
%  J_0to5=J.tot((1:50),1);
%  Vapp_0to5=Vapp(1,(1:50));
%  plot(Vapp_0to5,J_0to5);
%% Plot for electrons
%dfplot.JtotVapp(sol_CV_el, 0);
%set(gca,'YScale','log');

%% Electrostatic Plots for at starting, Vmax and Vmin
% dfplot.Vx(sol_CV,0);
% hold on
% dfplot.Vx(sol_CV,Vmaxtarray);
% hold on 
% dfplot.Vx(sol_CV,Vmintarray);
%% Anion Cation Plots for at Vmax 
%dfplot.acx(sol_CV,0);

% %%
% dfplot.acx(sol_CV,Vmaxtarray);
% %%
% dfplot.acx(sol_CV,Vmintarray);

%%
%dfplot.ELxnpxacx(sol_CV, 5)

%% To copy paste

% Scan rate = 0.0001 V/s
% Ionic species = 1
% 
% Scan rate = 100 V/s
% Ionic species = 2