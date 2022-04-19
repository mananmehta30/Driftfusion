function [J_electronic, delta_t, dV_by_dT_across_points, C_as_function_V_across_points] = capacitance_ana(sol_in)
par = sol_in.par;
[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_in);
J = dfana.calcJ(sol_in, "sub");%Get current density from dfana
J_tot=J.tot;
J_disp=J.disp;
J_electronic=J.tot-J.disp;
%% Get delta t
delta_t=t(:,2)-t(:,1); %change in time interval 
%% Create loop to calculate change in potential
for i=1:length(t)-1
    for j=1:length(x)
       dV_by_dT_across_points(i,j)= (V(i,j)-V(i+1,j))/delta_t;  %this is change in potential at each place for different times
    end
end
%%
k_scan = 0.001;
Vmax = 1.2;
Vmin = -1.2;
tpoints=(2*(Vmax-Vmin)/k_scan)+1;

time_array = [1,(Vmax/k_scan),3*(Vmax/k_scan)];
d=length(time_array);

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
legstr_c3 =[];
for i=1:length(time_array)
figure(122)
plot(x(par.pcum0(1,1):par.pcum0(1,3)),C_as_function_V_across_points(time_array(i),(par.pcum0(1,1):par.pcum0(1,3)))); 
legstr_c3{i} = ['Capacitance across insulator at t=', num2str(time_array(i))];
hold on
end
figure(122)
xlabel('Position [cm]')
ylabel('Capacitance (F/cm^2)')
legend(legstr_c3)
hold off
%% Capacitance vs voltage at the interface

V2=V;
V2(:,1)=[];%remove extra values
V2(481,:)=[];%remove extra values
%%
figure(123)
plot(V2(:,par.pcum0(3)),C_as_function_V_across_points(:,par.pcum0(3))); 
%legstr_ci3{i} = ['Capacitance across insulator at t=', num2str(time_array(i))];
%hold on
%end
figure(123)
xlabel('Voltage [V]')
ylabel('Capacitance (F/cm^2)')
%legend(legstr_ci3)
hold off
end
%% Get Debye Length
[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ions);
N_Debye=3 %calculate for three debye lengths
e = par.e;
V_T = par.kB*par.T;                     % Thermal votlage
epp_pvsk = e*par.epp0*par.epp(3);       % Perovskite absolute dielectric constant
N0 = par.Ncat(3);
debye_length = sqrt((epp_pvsk*V_T)/(e*N0));
%% Get the charge for the Debye length

total_electronic_charge=n+p;
total_ionic_charge=a+c;
electronic_charge_at_insulator_sc_interface =total_electronic_charge(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D));
ionic_charge_at_insulator_sc_interface =total_ionic_charge(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D));


%% Alternate way to calculate capacitance
% 1)Get parameters here %par = sol_in.par
% 2) extract various profiles[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_in);

%% Calculate charge

% 3) Way 1 to do it
%a) Integrate total charge Q_total=(p-n+Nd-Na) with change in potential
%across bulk and surface 
%fun = @(x) Q_total;
% E_s= sqrt( (-2*par.q)/par.epp*e0))integral(fun,V_bulk,V_surface)

