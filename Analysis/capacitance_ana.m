function [capacitance_device_electronic,capacitance_device_ionic] = capacitance_ana(sol_CV_with_ions)
  par_t = sol_CV_with_ions.par;
 [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ions);
%% Get Debye Length
% [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(soleq);
N_Debye=30; %calculate for specified debye lengths
e = par_t.e;
V_T = par_t.kB*par_t.T;                     % Thermal votlage
epp_pvsk = e*par_t.epp0*par_t.epp(3);       % Perovskite absolute dielectric constant
N0 = par_t.Ncat(3);
debye_length = sqrt((epp_pvsk*V_T)/(e*N0));
%% Get the charge for the Debye length
total_electronic_charge_density=(p-n)*e;%not n+p as they have different charges
total_ionic_charge_density=(c-a)*e;
x_perov_left = sol_CV_with_ions.par.dcum0(3);     % dcum is the device thickness
x_perov_right = sol_CV_with_ions.par.dcum0(4);
electronic_charge_at_insulator_sc_interface = total_electronic_charge_density(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
ionic_charge_at_insulator_sc_interface = total_ionic_charge_density(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
%% Intergrate to get space charge density from volumetric charge density
rho = dfana.calcrho(sol_CV_with_ions, "sub"); %Ask how to get the sub here
total_space_charge_per_unit_area=trapz(x, rho, 2);% Mistake is here, the charge integration gives one row for entire device
%% Find change in charge(s)
% For this you could use the MATLAB diff function
% https://uk.mathworks.com/help/matlab/ref/diff.html
del_q_ec=diff(electronic_charge_at_insulator_sc_interface);
del_q_ic=diff(ionic_charge_at_insulator_sc_interface);
del_sc=diff(total_space_charge_per_unit_area);
%% Find change in voltage applied
Vappt = dfana.calcVapp(sol_CV_with_ions);
for i=1:length(electronic_charge_at_insulator_sc_interface)-1
del_v(i)=Vappt(i+1)-Vappt(i);
end
%% Find change in voltage drop across the debye length
Vdrop=V(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
Vdrop2=diff(Vdrop);
% Vdrop3= V(:, x(sol_CV_with_ions.par.dcum0(3) sol_CV_with_ions.par.dcum0(3)+ N_Debye*debye_length))
Debye_left_V_index = find(x==sol_CV_with_ions.par.dcum0(3));
target=sol_CV_with_ions.par.dcum0(3)+ N_Debye*debye_length;
temp = abs(target - x);
closest = x(find(temp == min(abs(target - x))));
Debye_right_V_index=find(x==(closest));
Vdrop3=V(:,Debye_left_V_index)-V(:,Debye_right_V_index);
%potential across one side and then another
%plot the voltage and check the capacitance
%check same voltage at forward and reverse scan and compare (should be same as scan speed is slow)
%use dfplot.Qt to get electonic and ionic charge
%get the above values for the right hand side
%check if electronic and ionic capacitances add up to total capacitance
%freeze ions
%% Get total charge
% charge_debye = sum(total_space_charge_per_unit_area([Debye_left_V_index, Debye_right_V_index]),2);
charge_debye = sum(total_space_charge_per_unit_area(:,[Debye_left_V_index, Debye_right_V_index]),2)

%% Find C=change in charge by change in voltage
% for i=1:length(electronic_charge_at_insulator_sc_interface)-1
% capacitance_device_electronic(i)=(del_q_ec(i)/del_v(i))*e;
% capacitance_device_ionic(i)=(del_q_ic(i)/del_v(i))*e;
% end
%% Find capacitance across pvk layer

capacitance_device_electronic=(del_q_ec./Vdrop2);
capacitance_device_ionic=(del_q_ic./Vdrop2);
capacitance_del_sc=(del_sc./Vdrop2);

avg_c_electronic=mean(capacitance_device_electronic,2);
avg_c_ionic=mean(capacitance_device_ionic,2);
avg_c_spacecharge=mean(capacitance_del_sc,2);


avg_vdrop_electronic=mean(Vdrop2,2);
avg_vdrop_ionic=mean(Vdrop2,2);
%%
figure(1)
plot(Vdrop2, capacitance_device_electronic); 
xlabel('V drop')
ylabel('Electronic Capacitance at point across pvk layer with ions(F/cm^2)')

%%
figure(2)
plot(Vdrop2, capacitance_device_ionic); 
xlabel('V drop')
ylabel('Ionic Capacitance at point across pvk layer with ions(F/cm^2)')

%%
figure(3)
plot(avg_vdrop_electronic, avg_c_electronic); 
xlabel('Average Vdrop')
ylabel(' Average Electronic Capacitance across pvk layer  with ions(F/cm^2)')

%%
figure(4)
plot(avg_vdrop_ionic, avg_c_ionic); 
xlabel(' Average V drop')
ylabel('Average Ionic Capacitance across pvk layer with ions(F/cm^2)')

%%
figure(5)
plot(avg_vdrop_ionic,avg_c_spacecharge); 
xlabel(' Average V drop')
ylabel('Average space charge capacitance across pvk layer with ions(F/cm^2)')

%%
figure(6)
plot(avg_vdrop_ionic,del_sc); 
xlabel(' Average V drop')
ylabel('Average space charge across pvk layer with ions(Q/cm^2)')

%% Analytical solution Sze and Kwok

%% Define terms
% B = par.q/(par.k*par.t);
% n_p0=5;
% p_p0=5
% phi_s %Define surace potential
% F = sqrt((exp(-B*phi_s)+B*phi_s-1;

end
