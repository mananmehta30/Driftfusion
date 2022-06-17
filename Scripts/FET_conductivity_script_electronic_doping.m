%% Code pupose
% To get value of capacitance per area by integrating the current and using C(V)=J_displacement(V)/dV/dT
%Instead of using soleq.el better try with separate file

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');
par_for_ions = par_alox;     % Create temporary parameters object for overwriting parameters in loop
par_freeze_ions = par_alox;
%% Rough value of capacitance
A=1;
epsilon=par_for_ions.epp0*par_for_ions.epp(3)*par_for_ions.e;
d=par_for_ions.d(3);
Capacitance_rough=(A*epsilon)/d;
%% Set up parameters

Ncat_array = logspace(16, 19, 4);
kscan_index = [0.001;0.001;01;1];


%% Loop
for i = 1:length(Ncat_array)
    
    par_for_ions.Ncat(:) = Ncat_array(i);
    par_for_ions.Nani(:) = Ncat_array(i);
    disp(['Cation density = ', num2str(Ncat_array(i)), ' cm^-3']);
    par_freeze_ions.Ncat(:) = Ncat_array(i);
    par_freeze_ions.Nani(:) = Ncat_array(i);
    
    for j = 1:length(kscan_index) %loop to run for different electrode workfunction
        
              % This line is required to rebuild various arrays used DF
        par_for_ions.Phi_left = -4.9;
        disp(['LHS electrode workfunction = ', num2str(par_for_ions.Phi_left), ' eV']);
        par_for_ions.Phi_right = -4.9;
        disp(['RHS electrode workfunction = ', num2str(par_for_ions.Phi_right), ' eV']);

        
        par_freeze_ions.mu_c(:) = 0; %Freezing ions
        par_for_ions.mu_a(:) = 0;

        par_freeze_ions.Phi_left = -4.9;
        disp(['LHS electrode workfunction = ', num2str(par_for_ions.Phi_left), ' eV']);
        par_freeze_ions.Phi_right = -4.9;
        disp(['RHS electrode workfunction = ', num2str(par_for_ions.Phi_right), ' eV']);
        par_for_ions = refresh_device(par_for_ions);
        %% Find equilibrium
        soleq(i, j)= equilibrate(par_for_ions);
        soleq2(i, j)= equilibrate(par_freeze_ions);
        
 
        
        %% Current-voltage scan
        k_scan = kscan_index(j);
        disp(['Scan rate = ', num2str(kscan_index(j)), ' V/s']);
        Vmax = 5;
        Vmin = -5;
        tpoints=400;
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        
        sol_CV_with_ions(i, j) = doCV(soleq(i, j).ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
        sol_CV_without_ions(i, j) = doCV(soleq2(i, j).el, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
        %% Plot Vapp vs time
        % dfplot.Vappt(sol_CV)
        
        %% Plot JV scan
        %dfplot.JtotVapp(sol_CV, 0);
        %set(gca,'YScale','log')
        
        %% Plot anion and cation densities
        %dfplot.acx(sol_CV, 1/k_scan*[0:Vmax/3:Vmax]);
        
        %% Plot electron and hole profiles
        %dfplot.npx(sol_CV, 1/k_scan*[0:Vmax/3:Vmax]);
    end   
end




%%  Get and calculate Capacitance
Vappt = dfana.calcVapp(sol_CV_with_ions(i,j));
Vappt2 = dfana.calcVapp(sol_CV_without_ions(i,j));

[u,t,x,par_for_ions,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ions(i,j));
[u2,t2,x2,par_freeze_ions,dev2,n2,p2,a2,c2,V2] = dfana.splitsol(sol_CV_without_ions(i,j));

% Get Displacement current with and without ions
J_with_ions = dfana.calcJ(sol_CV_with_ions(i,j), "sub");
J_without_ions = dfana.calcJ(sol_CV_without_ions(i,j), "sub");

delta_Vapp=Vappt(:,2)-Vappt(:,1);%change in voltage dV
delta_t=t(:,2)-t(:,1);%change in time dt
delta_Vapp_by_delta_t=delta_Vapp/delta_t;%(dV/dt or k_scan basically)

J_disp_with_ions=J_with_ions.disp;
J_disp_with_ions=abs(J_disp_with_ions);%absolute values taken for clarity
C_with_ions=J_disp_with_ions/delta_Vapp_by_delta_t;

J_disp_without_ions=J_without_ions.disp;
J_disp_without_ions=abs(J_disp_without_ions);
C_without_ions=J_disp_without_ions/delta_Vapp_by_delta_t;
%% Plot Jdisp and Capacitance with ions across time

%Take central point of insulator
midpoint_insulator=(round((par_for_ions.pcum0(1)+par_for_ions.pcum0(2))/2));

figure(444)
plot(t,J_disp_with_ions(:,midpoint_insulator)); 
xlabel('Time')
ylabel('Displacement Current across insulator with ions')

figure(445)
plot(t,C_with_ions(:,midpoint_insulator)); 
xlabel('Time')
ylabel('Capacitance with ions')
%% Plot Jdisp and Capacitance without ions across time

%Take central point of insulator
midpoint_insulator=(round((par_for_ions.pcum0(1)+par_for_ions.pcum0(2))/2));
figure(544)
plot(t,J_disp_without_ions(:,midpoint_insulator)); 
xlabel('Time')
ylabel('Displacement Current across insulator without ions')

figure(545)
plot(t,C_without_ions(:,midpoint_insulator)); 
xlabel('Time')
ylabel('Capacitance without ions')
%% Capacitance_Manan
[CI, CWI] = capacitance_ana(sol_CV_with_ions);  
%% Call capacitance function
%[capacitance_device_electronic,capacitance_device_ionic] = capacitance_ana(sol_CV_with_ions);%call this function
[V, Q, C] = capacitance_ana_PC(sol_CV_with_ions, 2);    


%%

%%
% figure (499)
% plot(Vappt(2:end),capacitance_device_electronic); 
% xlabel('V applied')
% ylabel('Electronic Capacitance at point in an insulator with ions(F/cm^2)')
% %%
% figure (500)
% plot(Vappt(2:end), capacitance_device_ionic); 
% xlabel('V applied')
% ylabel('Ionic Capacitance at point in an insulator with ions(F/cm^2)')
%% Main Plots
% %% Plot Jdisp and Capacitance without ions across time
% 
% figure(446)
% plot(t,J_without_ions(:,122));
% xlabel('Time')
% ylabel('Displacement Current across insulator without ions')
% 
% figure(447)
% plot(t,C_without_ions(:,122));
% xlabel('Time')
% ylabel('Capacitance without ions')
% 
% %% Plot Capacitance vs Voltage applied for case with ions
% 
% figure(448)
% plot(Vappt,C_with_ions(:,122)); 
% xlabel('V applied')
% ylabel('Capacitance at point in an insulator with ions(F/cm^2)')
% 
% 
% %% Plot Capacitance vs Voltage applied for case without ions
% figure(449)
% plot(Vappt,C_without_ions(:,122)); 
% xlabel('V applied')
% ylabel('Capacitance at point in an insulator without ions(F/cm^2)')
% 
% %% Plot Capacitance vs Voltage applied for case with ions at interface
% 
% figure(450)
% plot(Vappt,C_with_ions(:,par.pcum0(1,3)+1)); 
% xlabel('V applied')
% ylabel('Capacitance at point in an insulator with ions at interface(F/cm^2)')
% 
% 
% %% Plot Capacitance vs Voltage applied for case without ions at interface
% figure(451)
% plot(Vappt,C_without_ions(:,par.pcum0(1,3)+1)); 
% xlabel('V applied')
% ylabel('Capacitance at point in an insulator without ions at interface(F/cm^2)')
% 
% %% df plot loop
% 
% for i=1:3
%     dfplot.ELnpx(sol_CV_without_ions,i+419)
%     hold on
% end
% hold off



%%
%Alternate ways to calculate capacitance

%1) https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=923259
%get surface electric field [FV, Frho] = dfana.calcF(sol_CV, "whole")(not
%working?)

