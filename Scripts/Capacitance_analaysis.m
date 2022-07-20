%% Code pupose and issue
% To get capacitance for different scan rates and other parameters.
% Code should store the Capacitance values for different voltage points in
% a array that can be accessed via data structures that store total,ionic
% and electronic capacitances

%Issue: The electronic capacitances seem to have some mini spikes and bumps
%that I am unable to understand as such (would changing tpoints with each loop help?)
% Also wanted to confirm if the data structures have been assigned properly.
% In capacitance_ana_PC some dot indexing issue is occuring

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');
par_for_ions = par_alox;     % Create temporary parameters object for overwriting parameters in loop
par_freeze_ions = par_alox;
%% Rough value of capacitance for MAPI, insulator and total
A=1;
epsilon_MAPI=par_for_ions.epp0*par_for_ions.epp(3)*par_for_ions.e;
d_MAPI=par_for_ions.d(3);
Capacitance_rough_MAPI=(A*epsilon_MAPI)/d_MAPI;
A=1;
epsilon_insulator=par_for_ions.epp0*par_for_ions.epp(2)*par_for_ions.e;
d_insulator=par_for_ions.d(1);
Capacitance_rough_insulator=(A*epsilon_insulator)/d_insulator;
Capacitance_total=1./((1/Capacitance_rough_MAPI)+(1/Capacitance_rough_insulator));
%% Set up parameters

Ncat_array = logspace(16, 19, 4);
kscan_array = [0.001;0.01;0.1];


%% Loop
for i = 1:length(Ncat_array)
    
    par_for_ions.Ncat(:) = Ncat_array(i);
    par_for_ions.Nani(:) = Ncat_array(i);
    disp(['Cation density = ', num2str(Ncat_array(i)), ' cm^-3']);
%     par_freeze_ions.Ncat(:) = Ncat_array(i);
%     par_freeze_ions.Nani(:) = Ncat_array(i);
    
    for j = 1:length(kscan_array) %loop to run for different k_scans
        
             
        par_for_ions.Phi_left = -4.9;
        disp(['LHS electrode workfunction = ', num2str(par_for_ions.Phi_left), ' eV']);
        par_for_ions.Phi_right = -4.9;
        disp(['RHS electrode workfunction = ', num2str(par_for_ions.Phi_right), ' eV']);

        
       % par_freeze_ions.mu_c(:) = 0; %Freezing ions
        %par_for_ions.mu_a(:) = 0;
        %par_freeze_ions.Phi_left = -4.9;
%   disp(['LHS electrode workfunction = ', num2str(par_for_ions.Phi_left), ' eV']);
%        par_freeze_ions.Phi_right = -4.9;
%        disp(['RHS electrode workfunction = ', num2str(par_for_ions.Phi_right), ' eV']);
       
par_for_ions = refresh_device(par_for_ions); % This line is required to rebuild various arrays used DF
        %% Find equilibrium
        soleq(i, j)= equilibrate(par_for_ions);
        %soleq2(i, j)= equilibrate(par_freeze_ions);
        
 
        
        %% Current-voltage scan
        k_scan = kscan_array(j);
        disp(['Scan rate = ', num2str(kscan_array(j)), ' V/s']);
        Vmax = 5;
        Vmin = -5;
        tpoints=400; 
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        
        sol_CV_with_ions(i, j) = doCV(soleq(i, j).ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
       % sol_CV_without_ions(i, j) = doCV(soleq2(i, j).el, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
      
    end   
end


%% Preallocate structures
for i=1:length(Ncat_array)
       C_ionic(i).N_cat=struct;
       C_electronic(i).N_cat=struct;
       C_total(i).N_cat=struct;
     
end
%% Preallocate structures to store k_scan values

for i=1:length(Ncat_array)
    
 for j=1:length(kscan_array)
       C_ionic(i).N_cat(j).k_scan=struct;
       C_electronic(i).N_cat(j).k_scan=struct;
       C_total(i).N_cat(j).k_scan=struct;
 end
 end


%% Capacitance_Manan Analysis

Vappt = dfana.calcVapp(sol_CV_with_ions(1,1)); 


for i = 1:length(Ncat_array)
    for j = 1:length(kscan_array)
       
        [Ctotal,Celectronic,Cionic] = capacitance_ana(sol_CV_with_ions(i,j),Vappt);%call this function
       
        C_ionic(i).N_cat(j).k_scan=Cionic;
        C_electronic(i).N_cat(j).k_scan=Celectronic;
        C_total(i).N_cat(j).k_scan=Ctotal;
        
    end
end

 


%% Plot carrier concentration at interface as function Vapp for different ion densities
scanrate_index = 1;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    figure(203)
    plot(Vappt, C_total(i).N_cat(scanrate_index).k_scan)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i)), 'cm-3'];
    hold on
end


figure(203)
xlabel('Voltage [V]')
ylabel('Total Capacitance [cm-3]')
legend(legstr_n3)
hold off

for i = 1:length(Ncat_array)
    figure(204)
    plot(Vappt, C_ionic(i).N_cat(scanrate_index).k_scan)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i)), 'cm-3'];
    hold on
end


figure(204)
xlabel('Voltage [V]')
ylabel('Ionic Capacitance [cm-3]')
legend(legstr_n3)
hold off

for i = 1:length(Ncat_array)
    figure(205)
    plot(Vappt, C_electronic(i).N_cat(scanrate_index).k_scan)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i)), 'cm-3'];
    hold on
end


figure(205)
xlabel('Voltage [V]')
ylabel('Electronic Capacitance [cm-3]')
legend(legstr_n3)
hold off

%% Plot carrier concentration at interface as function Vapp for different scan rates
ion_concentration_index = 1;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(kscan_array)
    figure(203)
    plot(Vappt, C_total(ion_concentration_index).N_cat(i).k_scan)
    legstr_n3{i} = ['kscan=', num2str(kscan_array(i)), 'V/s'];
    hold on
end


figure(203)
xlabel('Voltage [V]')
ylabel('Total Capacitance [cm-3]')
legend(legstr_n3)
hold off

for i = 1:length(kscan_array)
    figure(204)
    plot(Vappt, C_ionic(ion_concentration_index).N_cat(i).k_scan)
    legstr_n3{i} = ['kscan=', num2str(kscan_array(i)), 'V/s'];
    hold on
end


figure(204)
xlabel('Voltage [V]')
ylabel('Ionic Capacitance [cm-3]')
legend(legstr_n3)
hold off

for i = 1:length(kscan_array)
    figure(205)
    plot(Vappt, C_electronic(ion_concentration_index).N_cat(i).k_scan)
    legstr_n3{i} = ['kscan=', num2str(kscan_array(i)), 'V/s'];
    hold on
end


figure(205)
xlabel('Voltage [V]')
ylabel('Electronic Capacitance [cm-3]')
legend(legstr_n3)
hold off
%% Call capacitance function 

%[V, Q, C] = capacitance_ana_PC(sol_CV_with_ions, 2);   
% figure(7445)
% plot(Vappt, abs(C_debye_layers), Vappt, abs(C_debye_electronic), '-.', Vappt, abs(C_debye_ionic), '--')
% legend('Total Capacitance','Electronic Capacitance','Ionic Capacitance')
% xlabel('Voltage applied')
% ylabel('Capacitances (F/cm^2)')




