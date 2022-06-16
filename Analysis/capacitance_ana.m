function [capacitance_device_electronic,capacitance_device_ionic] = capacitance_ana(sol_CV_with_ion)
par_t = sol_CV_with_ion.par;
[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ion);
%% Get Debye Length
% [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(soleq);
N_Debye=2; %calculate for specified debye lengths
e = par_t.e;
V_T = par_t.kB*par_t.T;                     % Thermal votlage
epp_pvsk = e*par_t.epp0*par_t.epp(3);       % Perovskite absolute dielectric constant
N0 = par_t.Ncat(3);
debye_length = sqrt((epp_pvsk*V_T)/(e*N0));
%% Get the charge for the Debye length
rho = dfana.calcrho(sol_CV_with_ion, "whole");
total_space_charge_per_unit_area = e.*trapz(x, rho, 2);
total_electronic_charge_density=(p-n)*e;%not n+p as they have different charges
total_ionic_charge_density=(c-a)*e;
% PC - these calculations of charge miss out the static background dopants
% and counter-ionic charges for the mobile ionic species. Better to use
% dfana.calcrho(sol) - I see you have done this below in fact

x_perov_left = par_t.dcum0(3);     % dcum is the device thickness
x_perov_right = par_t.dcum0(4);

electronic_charge_at_insulator_sc_interface = total_electronic_charge_density(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
ionic_charge_at_insulator_sc_interface = total_ionic_charge_density(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
%% Find change in charge(s)
% For this you could use the MATLAB diff function
% https://uk.mathworks.com/help/matlab/ref/diff.html
del_q_ec=diff(electronic_charge_at_insulator_sc_interface);
del_q_ic=diff(ionic_charge_at_insulator_sc_interface);
del_sc=diff(total_space_charge_per_unit_area);
%% Find change in voltage applied
Vappt = dfana.calcVapp(sol_CV_with_ion);
del_v = diff(Vappt);
% PC - the built-in function DIFF can be used here to make it more concise-
% I don't think you need the differences however as you can use the
% gradient function instead see below

%% Find change in voltage drop across the debye length(two methods
Vdrop = V(:, x > x_perov_left & x < x_perov_left + N_Debye*debye_length);
Vdrop2 = diff(Vdrop);
Vdrop = V(:, dsearchn(x', x_perov_left)) - V(:, dsearchn(x', x_perov_left + N_Debye*debye_length));

%%%% Method 1 to get difference in voltage
Debye_left_V_index = find(x==sol_CV_with_ion.par.dcum0(3));
target=sol_CV_with_ion.par.dcum0(3)+ N_Debye*debye_length;
temp = abs(target - x);
closest = x(find(temp == min(abs(target - x))));
% PC- blimey- this looks complicated! I think you can use DSEARCHN(X, X_REQUESTED) to get
% the nearest point index to X_REQUESTED

Debye_right_V_index=find(x==(closest));
Vdrop3=V(:,Debye_left_V_index)-V(:,Debye_right_V_index);

%%%% Method 3 to get difference in voltage
deltaV = dfana.deltaVt(sol_CV_with_ion, Debye_left_V_index, Debye_right_V_index);

%potential across one side and then another
%plot the voltage and check the capacitance
%check same voltage at forward and reverse scan and compare (should be same as scan speed is slow)
%use dfplot.Qt to get electonic and ionic charge
%get the above values for the right hand side
%check if electronic and ionic capacitances add up to total capacitance
%freeze ions

%% Get rho across debye length
pl = Debye_left_V_index;
pr = Debye_right_V_index;

rho = dfana.calcrho(sol_CV_with_ion,"whole");
% n0=n(:, pl:pr);
% p0=p(:, pl:pr);
% a0=a(:, pl:pr);
% c0=c(:, pl:pr);
NA0 = repmat(dev.NA, length(t), 1);
ND0 = repmat(dev.ND, length(t), 1);
Nani0 = repmat(dev.Nani, length(t), 1);
Ncat0 = repmat(dev.Ncat, length(t), 1);

rho0 = rho(:, pl:pr); % + par.z_a*a + par.z_c*c - par.z_c*Nani0 - par.z_c*Ncat0;

% PC- was RHO0 a check?
total_space_charge_per_unit_area0=e.*trapz(x(pl:pr), rho(:, pl:pr), 2);
% PC - I think I may have been too entusiastic when I told you to be
% explicit with variable names! e.g. this could be called SIGMA and then
% you could describe what it is in the comments

electronic_rho = -n + p - NA0 + ND0;
total_electronic_charge_per_unit_area0=e.*trapz(x(pl:pr), electronic_rho(:, pl:pr), 2);

ionic_rho = par.z_a*a + par.z_c*c - par.z_c*Nani0 - par.z_c*Ncat0;
total_ionic_charge_per_unit_area0=e.*trapz(x(pl:pr), ionic_rho(:, pl:pr), 2);

%% Space charge total
Q = e*trapz(x(pl:pr), rho(:, pl:pr), 2);

electronic_rho=p - n - NA0 + ND0;
ionic_rho= par.z_a*a + par.z_c*c - par.z_c*Nani0 - par.z_c*Ncat0;
% PC- you have already calculated the above a number of times?
% Again these should include the static dopant terms otherwise it doesn't
% make any sense- I have included these for you
EQ = e*trapz(x(pl:pr), electronic_rho(:, pl:pr), 2);
IQ= e*trapz(x(pl:pr), ionic_rho(:, pl:pr), 2);

%% Get capacitance as total charge (per cm2) divided by dV
%PC- remember C = dQdV in the most general form
% del_space_charge=diff(Q);
% del_V_debye=diff(deltaV);
C_debye_layers = gradient(Q, Vdrop);

%del_space_electronic=diff(EQ);
C_debye_electronic = gradient(EQ, Vdrop);

%del_ionic_charge=diff(IQ);
C_debye_ionic = gradient(IQ, Vdrop); %del_ionic_charge./del_V_debye;

% deltaV2=deltaV;
% deltaV2(1,:)=[];  % PC- when you use the GRADIENT function it also
% approximates the boundary values so that you don't lose points

%subplot(2,1,1);
figure(7444)
semilogy(Vdrop, abs(C_debye_layers), Vdrop, abs(C_debye_electronic), '-.', Vdrop, abs(C_debye_ionic), '--')
legend('Total Capacitance','Electronic Capacitance','Ionic Capacitance')
xlabel('V across debye layers')
ylabel('Capacitances (F/cm^2)')

%% Find C=change in charge by change in voltage
% for i=1:length(electronic_charge_at_insulator_sc_interface)-1
% capacitance_device_electronic(i)=(del_q_ec(i)/del_v(i))*e;
% capacitance_device_ionic(i)=(del_q_ic(i)/del_v(i))*e;
% end
figure(45)
Vappt2=Vappt;
% Vappt2(:,1)=[];
Vappt3=transpose(Vappt2);
plot(Vappt3,C_debye_layers)
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
plot(Vdrop2, abs(capacitance_device_electronic)); 
xlabel('V drop')
ylabel('Electronic Capacitance at point across pvk layer with ions(F/cm^2)')

%%
figure(2)
plot(Vdrop2, abs(capacitance_device_ionic)); 
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
%%
% figure(7)
% plot(Vdrop3,total_electronic_charge_density_inthedebye); 
% xlabel('V drop across Debye')
% ylabel('Total volumetric charge density across pvk layer with ions(Q/cm^2)')

%% Analytical solution Sze and Kwok

%% Define terms
% B = par.q/(par.k*par.t);
% n_p0=5;
% p_p0=5
% phi_s %Define surace potential
% F = sqrt((exp(-B*phi_s)+B*phi_s-1;

end
