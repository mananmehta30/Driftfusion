%% Code pupose
% To get value of capacitance per area by integrating the current and using C(V)=J_displacement(V)/dV/dT
%Instead of using soleq.el better try with separate file

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');
par = par_alox;     % Create temporary parameters object for overwriting parameters in loop
%% Rough value of capacitance
A=1;
epsilon=par.epp0*par.epp(3)*par.e;
d=par.d(3);
Capacitance_rough=(A*epsilon)/d;
%% Par file for ion case
par.Ncat(:) = 1e19; %Simulating for commonly reported ionic density
par.Nani(:) = 1e19; %Simulating for commonly reported ionic density 
disp(['Cation density = ', num2str(par.Ncat(3)), ' cm^-3']);%num2str=Convert numbers to character representation
      
par.Phi_left = -4.9;
disp(['LHS electrode workfunction = ', num2str(par.Phi_left), ' eV']);
par.Phi_right = -4.9;
disp(['RHS electrode workfunction = ', num2str(par.Phi_right), ' eV']);    


%% Par file for no ion case

par2 = par_alox;
par2.mu_c(:) = 0; %Freezing ions
par.mu_a(:) = 0;

par2.Phi_left = -4.9;
disp(['LHS electrode workfunction = ', num2str(par.Phi_left), ' eV']);
par2.Phi_right = -4.9;
disp(['RHS electrode workfunction = ', num2str(par.Phi_right), ' eV']);

%% Find equilibrium
soleq= equilibrate(par);
soleq2= equilibrate(par2);
%% Current-voltage scan
k_scan = 0.001;
Vmax = 1.2;
Vmin = -Vmax;
tpoints=241;
%Try freezing the ions (keep mobility =0) and then run sol_CV to check
%later on
%sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_with_ions = doCV(soleq.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
sol_CV_without_ions = doCV(soleq2.el, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);

%%  Get and calculate Capacitance
Vappt = dfana.calcVapp(sol_CV_with_ions);
Vappt2 = dfana.calcVapp(sol_CV_without_ions);

[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ions);
[u2,t2,x2,par2,dev2,n2,p2,a2,c2,V2] = dfana.splitsol(sol_CV_without_ions);

% Get Displacement current with and without ions
J_with_ions = dfana.calcJ(sol_CV_with_ions, "sub");
J_without_ions = dfana.calcJ(sol_CV_without_ions, "sub");

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
midpoint_insulator=(round((par.pcum0(1)+par.pcum0(2))/2));

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
midpoint_insulator=(round((par.pcum0(1)+par.pcum0(2))/2));
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

