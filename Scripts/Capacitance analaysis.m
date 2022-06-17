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
kscan_index = [0.01;0.1;1];


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
      
    end   
end




%%  Get and calculate Capacitance
% Vappt = dfana.calcVapp(sol_CV_with_ions(i,j));
% Vappt2 = dfana.calcVapp(sol_CV_without_ions(i,j));
% 
% [u,t,x,par_for_ions,dev,n,p,a,c,Ctotal] = dfana.splitsol(sol_CV_with_ions(i,j));
% [u2,t2,x2,par_freeze_ions,dev2,n2,p2,a2,c2,V2] = dfana.splitsol(sol_CV_without_ions(i,j));

% % Get Displacement current with and without ions
% J_with_ions = dfana.calcJ(sol_CV_with_ions(i,j), "sub");
% J_without_ions = dfana.calcJ(sol_CV_without_ions(i,j), "sub");
% 
% delta_Vapp=Vappt(:,2)-Vappt(:,1);%change in voltage dV
% delta_t=t(:,2)-t(:,1);%change in time dt
% delta_Vapp_by_delta_t=delta_Vapp/delta_t;%(dV/dt or k_scan basically)
% 
% J_disp_with_ions=J_with_ions.disp;
% J_disp_with_ions=abs(J_disp_with_ions);%absolute values taken for clarity
% C_with_ions=J_disp_with_ions/delta_Vapp_by_delta_t;
% 
% J_disp_without_ions=J_without_ions.disp;
% J_disp_without_ions=abs(J_disp_without_ions);
% C_without_ions=J_disp_without_ions/delta_Vapp_by_delta_t;
% %% Plot Jdisp and Capacitance with ions across time
% 
% %Take central point of insulator
% midpoint_insulator=(round((par_for_ions.pcum0(1)+par_for_ions.pcum0(2))/2));
% 
% figure(444)
% plot(t,J_disp_with_ions(:,midpoint_insulator)); 
% xlabel('Time')
% ylabel('Displacement Current across insulator with ions')
% 
% figure(445)
% plot(t,C_with_ions(:,midpoint_insulator)); 
% xlabel('Time')
% ylabel('Capacitance with ions')
% %% Plot Jdisp and Capacitance without ions across time
% 
% %Take central point of insulator
% midpoint_insulator=(round((par_for_ions.pcum0(1)+par_for_ions.pcum0(2))/2));
% figure(544)
% plot(t,J_disp_without_ions(:,midpoint_insulator)); 
% xlabel('Time')
% ylabel('Displacement Current across insulator without ions')
% 
% figure(545)
% plot(t,C_without_ions(:,midpoint_insulator)); 
% xlabel('Time')
% ylabel('Capacitance without ions')
%% Create data structure to store capacitance value
C = struct;
C.C_total=zeros(length(Ncat_array),length(kscan_index));
C.C_ionic=zeros(length(Ncat_array),length(kscan_index));
C.C_electronic=zeros(length(Ncat_array),length(kscan_index));

%% Capacitance_Manan Analysis
Vappt = dfana.calcVapp(sol_CV_with_ions(1,1)); 


for i = 1:length(Ncat_array)
    for j = 1:length(kscan_index)
       
        [Ctotal, Cionic, Celectronic] = capacitance_ana(sol_CV_with_ions(i,j),Vappt);%call this function
        C.C_total(i,j)=Ctotal;
        C.C_ionic(i,j)=Cionic;
        C.C_ionic(i,j)=Celectronic;
    end
end

%% Call capacitance function
%[capacitance_device_electronic,capacitance_device_ionic] = capacitance_ana(sol_CV_with_ions);%call this function
[V, Q, C] = capacitance_ana_PC(sol_CV_with_ions, 2);    


%% Plots

figure(7445)
plot(Vappt, abs(C_debye_layers), Vappt, abs(C_debye_electronic), '-.', Vappt, abs(C_debye_ionic), '--')
legend('Total Capacitance','Electronic Capacitance','Ionic Capacitance')
xlabel('Voltage applied')
ylabel('Capacitances (F/cm^2)')





%%
%Alternate ways to calculate capacitance

%1) https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=923259
%get surface electric field [FV, Frho] = dfana.calcF(sol_CV, "whole")(not
%working?)

