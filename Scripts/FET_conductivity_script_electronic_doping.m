

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
soleq= equilibrate(par);
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
%% Get delta t
% delta_t=t(:,2)-t(:,1); %change in time interval 
%% Create loop to calculate change in potential
for i=1:length(t)-1
    for j=1:length(x)
        [J_electronic, delta_t, dV_by_dT_across_points, C_as_function_V_across_points] = capacitance_ana(sol_CV(i,j));  %this is change in potential at each place for different times
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
J_disp(4801,:)=[]; 
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
%Alternate ways to calculate capacitance

%1) https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=923259
%get surface electric field dfplot.ft

dfplot.Ft(sol_CV,3);