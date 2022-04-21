function [capacitance_device_electronic,capacitance_device_ionic] = capacitance_ana(sol_CV_with_ions)
  par_t = sol_CV_with_ions.par;
 [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ions);
%% Get Debye Length
% [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(soleq);
N_Debye=2; %calculate for specified debye lengths
e = par_t.e;
V_T = par_t.kB*par_t.T;                     % Thermal votlage
epp_pvsk = e*par_t.epp0*par_t.epp(3);       % Perovskite absolute dielectric constant
N0 = par_t.Ncat(3);
debye_length = sqrt((epp_pvsk*V_T)/(e*N0));
%% Get the charge for the Debye length
total_electronic_charge=n+p;
total_ionic_charge=a+c;
x_perov_left = sol_CV_with_ions.par.dcum0(3);     % dcum is the device thickness
x_perov_right = sol_CV_with_ions.par.dcum0(4);
electronic_charge_at_insulator_sc_interface = total_electronic_charge(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
ionic_charge_at_insulator_sc_interface = total_ionic_charge(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
%% Intergrate to get space charge density from volumetric charge density
rho = dfana.calcrho(sol_CV_with_ions, "whole");
elec_space_charge=trapz(x, rho, 2);
%% Find change in charge(s)
% For this you could use the MATLAB diff function
% https://uk.mathworks.com/help/matlab/ref/diff.html
del_q_ec=diff(electronic_charge_at_insulator_sc_interface);
del_q_ic=diff(ionic_charge_at_insulator_sc_interface);
del_sc=diff(elec_space_charge);
%% Find change in voltage applied
Vappt = dfana.calcVapp(sol_CV_with_ions);
for i=1:length(electronic_charge_at_insulator_sc_interface)-1
del_v(i)=Vappt(i+1)-Vappt(i);
end
%% Find change in voltage drop across the debye length
Vdrop=V(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
Vdrop2=diff(Vdrop);
%% Find C=change in charge by change in voltage
% for i=1:length(electronic_charge_at_insulator_sc_interface)-1
% capacitance_device_electronic(i)=(del_q_ec(i)/del_v(i))*e;
% capacitance_device_ionic(i)=(del_q_ic(i)/del_v(i))*e;
% end
%% Find capacitance across pvk layer

capacitance_device_electronic=(del_q_ec./Vdrop2)*e;
capacitance_device_ionic=(del_q_ic./Vdrop2)*e;


%%
figure(1)
plot(Vdrop2, capacitance_device_electronic); 
xlabel('Vdrop')
ylabel('Electronic Capacitance at point in an insulator with ions(F/cm^2)')

%%
figure(2)
plot(Vdrop2, capacitance_device_ionic); 
xlabel('V drop')
ylabel('Ionic Capacitance at point in an insulator with ions(F/cm^2)')

%% Analytical solution Sze and Kwok

%% Define terms
% B = par.q/(par.k*par.t);
% n_p0=5;
% p_p0=5
% phi_s %Define surace potential
% F = sqrt((exp(-B*phi_s)+B*phi_s-1;

end
