%% Code pupose
% To get value of capacitance per area by integrating the current and using C(V)=J_displacement(V)/dV/dT
%Instead of using soleq.el better try with separate file

%% Initialize driftfusion
initialise_df
%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');

par = par_alox;     % Create temporary parameters object for overwriting parameters in loop
par.Ncat(:) = 1e19; %Simulating for commonly reported ionic density
par.Nani(:) = 1e19; %Simulating for commonly reported ionic density 

par2 = par_alox;
par2.N_ionic_species=0; %keeping no ionic density
par2.Ncat(:) = 0; %Simulating for 0 ionic density
par2.Nani(:) = 0; %Simulating for 0 ionic density 

disp(['Cation density = ', num2str(par.Ncat(3)), ' cm^-3']);%num2str=Convert numbers to character representation
      
par.Phi_left = -4.9;
disp(['LHS electrode workfunction = ', num2str(par.Phi_left), ' eV']);
par.Phi_right = -4.9;
disp(['RHS electrode workfunction = ', num2str(par.Phi_right), ' eV']);    

par2.Phi_left = -4.9;
disp(['LHS electrode workfunction = ', num2str(par.Phi_left), ' eV']);
par2.Phi_right = -4.9;
disp(['RHS electrode workfunction = ', num2str(par.Phi_right), ' eV']);

%% Find equilibrium
soleq= equilibrate(par);
soleq2= equilibrate(par2);
%% Current-voltage scan
k_scan = 0.001;
Vmax = 10;
Vmin = -Vmax;
tpoints=(2*(Vmax-Vmin)/(10*k_scan))+1;

%sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_with_ions = doCV(soleq.ion, 0, 0, Vmax, Vmin, k_scan, 1, 1000);
sol_CV_without_ions = doCV(soleq2.el, 0, 0, Vmax, Vmin, k_scan, 1, 1000);

%% Vappt and other parameters
Vappt = dfana.calcVapp(sol_CV_with_ions);
Vappt2 = dfana.calcVapp(sol_CV_without_ions);

[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ions);
[u2,t2,x2,par2,dev2,n2,p2,a2,c2,V2] = dfana.splitsol(sol_CV_without_ions);

%% Get J disp with and without ions
% J_with_ions = dfana.calcJ(sol_CV_with_ions, "sub");
% J_without_ions = dfana.calcJ(sol_CV_without_ions, "sub");
% 
% delta_Vapp=Vappt(:,2)-Vappt(:,1);%change in voltage dV
% delta_t=t(:,2)-t(:,1);%change in time dt
% delta_Vapp_by_delta_t=delta_Vapp/delta_t;%(dV/dt or k_scan basically)
% 
% J_with_ions=J_with_ions.disp;
% J_with_ions=abs(J_with_ions);%absolute values taken for clarity
% C_with_ions=J_with_ions/delta_Vapp_by_delta_t;
% 
% J_without_ions=J_without_ions.disp;
% J_without_ions=abs(J_without_ions);
% C_without_ions=J_without_ions/delta_Vapp_by_delta_t;
% %% Plot Jdisp and Capacitance with ions across time
% figure(444)
% plot(t,J_with_ions(:,122)); 
% xlabel('Time')
% ylabel('Displacement Current across insulator with ions')
% 
% figure(445)
% plot(t,C_with_ions(:,122)); 
% xlabel('Time')
% ylabel('Capacitance with ions')
%% Call capacitance function
[capacitance_device_electronic,capacitance_device_ionic] = capacitance_ana(sol_CV_with_ions);%call this function
       
%%

%%
plot(Vappt2,capacitance_device_electronic); 
xlabel('V applied')
ylabel('ElectronicCapacitance at point in an insulator with ions(F/cm^2)')
%%
plot(Vappt2,capacitance_device_ionic); 
xlabel('V applied')
ylabel('Ionic Capacitance at point in an insulator with ions(F/cm^2)')
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

