

%% Code pupose
% To get value of capacitance per area by integrating the current and using
% C(V)=J(V)/dV/dT

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');

par = par_alox;     % Create temporary parameters object for overwriting parameters in loop


par.Ncat(:) = 1e19; %Simulating for commonly reported ionic density
par.Nani(:) = 1e19; %Simulating for commonly reported ionic density 
   
disp(['Cation density = ', num2str(par.Ncat(3)), ' cm^-3']);%num2str=Convert numbers to character representation
   
        
par.Phi_left = -4.9;
disp(['LHS electrode workfunction = ', num2str(par.Phi_left), ' eV']);
par.Phi_right = -4.9;
disp(['RHS electrode workfunction = ', num2str(par.Phi_right), ' eV']);     

%% Find equilibrium
soleq= equilibrate(par);% Find and 
%% Current-voltage scan
k_scan = 0.001;
Vmax = 1.2;
Vmin = -1.2;
tpoints=(2*(Vmax-Vmin)/k_scan)+1;

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
  sol_CV = doCV(soleq.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);

        %% Plot Vapp vs time
         %dfplot.Vappt(sol_CV)
        
        %% Plot JV scan
        dfplot.JtotVapp(sol_CV, 0);
        set(gca,'YScale','log')
        
        %% Plot anion and cation densities
        dfplot.acx(sol_CV, 1/k_scan*[0:Vmax:3*Vmax]);
        
        %% Plot electron and hole profiles
        dfplot.npx(sol_CV, (1/k_scan)*Vmax);

        %% Plot current as a function of time
        dfplot.Jt(sol_CV,100);
%% Capacitance calculation
% To get value of capacitance per area by integrating the current and using
% C(V)=J(V)/dV/dT
[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV);
J = dfana.calcJ(sol_CV, "sub");%Get current density from dfana
J_tot=J.tot;
J_disp=J.disp;
J_electronic=J.tot-J.disp;
%% Get delta t
delta_t=t(:,2)-t(:,1); %change in time interval 
%% Create loop to calculate change in potential
for i=1:length(t)-1
    for j=1:length(x)
        [capacitance] = capacitance_ana(sol_CV(i,j))  %this is change in potential at each place for different times
    end
end
% dV_by_dT_across_points(i,j) = delta_v_across_points(i,j)/delta_t; %dividing it by the time interval to get dV/dt
% %% Create loop to calculate change in potential
% for i=1:length(t)-1
%     for j=1:length(x)
%         delta_v_across_points(i,j)= V(i,j)-V(i+1,j);  %this is change in potential at each place for different times
%     end
% end
% dV_by_dT_across_points(i,j) = delta_v_across_points(i,j)/delta_t; %dividing it by the time interval to get dV/dt
% 
% %%

%% Remove the first column since it does not contribute
dV_by_dT_across_points(:,1) = []; 
%% Remove final row since it does not contribute
J_disp(481,:)=[]; 
%% Calculate capacitance at each point
for i=1:length(t)-1
    for j=1:length(x)-1
        C_as_function_V_across_points(i,j) = J_disp(i,j)/dV_by_dT_across_points(i,j);%Get C=J(V)/(dV/dt)
    end
end

%% Plot capacitance as a function of position
x(:,1)=[]; %Remove the first point to since it doesnt come inside the calculation
figure(7464)
plot(x(100:320),C_as_function_V_across_points(6,(100:320))); 
%%
for i = 1:length(Ncat_array)
    figure(100)
    semilogy(workfunction_LHS, sigma_n_bar_VpeakM(i, :))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
    xlim([workfunction_LHS(1) workfunction_LHS(length(workfunction_LHS))])
    
end  

for i = 1:length(Ncat_array)
    figure(101)
    semilogy(workfunction_LHS, sigma_p_bar_VpeakM(i, :))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak hole conductivity [S cm-1]')
    legstr_p{i} = ['Ncat =', num2str(Ncat_array(i))];
    xlim([workfunction_LHS(1) workfunction_LHS(length(workfunction_LHS))])
    
end  
figure(100)
legend(legstr_n)
hold off
figure(101)
legend(legstr_p)
hold off

%% Plot cation density as a function of voltage at the interface
cation_index=3;
legstr_n3 =[];
for i = 1:length(workfunction_LHS)
    cation_density_interface = sol_CV(cation_index, i).u(:, par.pcum0(3)+1,4); %Whats the value here for cation? Should par.pcum0(3) be different?
    figure(205)
    semilogy(Vappt, cation_density_interface)
    legstr_n3{i} = ['Phi =', num2str(workfunction_LHS(i))];
    hold on
end

figure(205)
xlabel('Voltage [V]')
ylabel('Cation density interface (cm-3)')
legend(legstr_n3)
hold off
%% Plot average conductivity
for j = 1:length(workfunction_LHS)
    figure(201)
    semilogy(Vappt, squeeze(sigma_n_barM(3, j, :))) %remove singleton dimensions (reason?)
    legstr_n2{j} = ['\Phi_l =', num2str(workfunction_LHS(j))];
    xlim([Vmin Vmax]);
    hold on
end

for j = 1:length(workfunction_LHS)
    figure(202)
    semilogy(Vappt, squeeze(sigma_p_barM(3, j, :)))
    legstr_p2{j} = ['\Phi_l =', num2str(workfunction_LHS(j))];
    xlim([Vmin Vmax]);
    hold on
end

figure(201)
xlabel('Voltage [V]')
ylabel('Average electron conductivity [Siemens]')
legend(legstr_n2)
hold off

figure(202)
xlabel('Voltage [V]')
ylabel('Average hole conductivity [Siemens]')
legend(legstr_p2)
hold off


%% Plot carrier concentration at interface as function Vapp for different ion densities
workfunction_index =1;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    n_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3), 2); %n_int takes the electron density at the interface for all times, at space
    figure(203)
    semilogy(Vappt, n_int)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

for i = 1:length(Ncat_array)
    p_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3), 3);%par.pcum0(3) that corresponds to MAPI for variable number 2 that is the hole density. 

    figure(204)
    semilogy(Vappt, p_int)
    legstr_p3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

figure(203)
xlabel('Voltage [V]')
ylabel('electron density interface (cm-3)')
legend(legstr_n3)
hold off

figure(204)
xlabel('Voltage [V]')
ylabel('hole density interface (cm-3)')
legend(legstr_p3)
hold off

%% Plot electron and hole profiles at Vmax as a function of position
workfunction_index =1;
legstr_npx = {'', '', ''};
Vmax_t_index = find(sol_CV(1, 1).t == 80);
V_INDEX=Vmax/k_scan; %Better way to use the index?
for i = 1:length(Ncat_array)
    dfplot.npx(sol_CV(i, workfunction_index), Vmax/k_scan);% Vmax/k_scan) %this is the time index if I am correct
    legstr_npx{2*i-1 + 3} = ['n, Ncat =', num2str(Ncat_array(i))];
    legstr_npx{2*i + 3} = ['p, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_npx)
ylim([1e-1, 1e12])

%% Plot cation and anion densities at Vmax as a function of position
legstr_acx = {'', '', ''};
for i = 1:length(Ncat_array)
    dfplot.acx(sol_CV(i, workfunction_index), Vmax/k_scan)
    legstr_acx{2*i-1 + 3} = ['a, Ncat =', num2str(Ncat_array(i))];
    legstr_acx{2*i + 3} = ['c, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_acx)
%ylim([1e-1, 1e12])

%% Plot potential as a function position
workfunction_index =1;
legstr_Vx = {'dielectric', 'interface', 'perovskite'};
for i = 1:length(Ncat_array)
    dfplot.Vx(sol_CV(i, workfunction_index), 0);%Vmax/k_scan)
    legstr_Vx{i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_Vx)
%ylim([1e-1, 1e12])

%% Trying some dfplots
%   
for i = 4:-1:1
    dfplot.Vxacx(sol_CV(i,1), 0)
    subplot(2,1,1);
    hold on
    subplot(2,1,2)
    hold on
end
 subplot(2,1,1);
    hold off
    subplot(2,1,2)
    hold off
%Why does the graph disappear?

%% Potential due to ionic effect

   % dfplot.Vionxacx(sol_CV(1,7), 0)
   
   %% DFPLOT
  
   
dfplot.acx(sol_CV(3, 1),180);
%this is not time. check again
%% Plot individual values
%Ask how makemovie works
%  makemovie(sol_CV, @dfplot.npx, 0, [0, 1.5e18], 'npx', true, true);


%% Conductivity profiles
% So systematically you could look at the following.

% 1)Electrode workfunctions
% 2)Ion density 
% 3)1 ion, opposite charge (i.e. mobile anions)
% 4)2 ions
% 5)Different ion densities

% Ideally you would set some of these parameter explorations up as loops and extract peak conductivity then plot 
% on a contour plot with x = Ion density, y = Electrode workfunctions, z = peak conductivity for example.
%% Questions

%Why is dfana used for generating the voltage steps? Can it be done
%directly instead?
