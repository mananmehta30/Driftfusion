function [C_debye_layers,C_debye_electronic,C_debye_ionic] = capacitance_ana(sol_CV_with_ions,Vappt)
par_t = sol_CV_with_ions.par;
[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV_with_ions);
%% Get Debye Length
N_Debye=5; %calculate for specified debye lengths
e = par_t.e;
V_T = par_t.kB*par_t.T;                     % Thermal votlage
epp_pvsk = e*par_t.epp0*par_t.epp(3);       % Perovskite absolute dielectric constant
N0 = par_t.Ncat(3);
debye_length = sqrt((epp_pvsk*V_T)/(e*N0));
%% Get beginning poinnt of the Perovskite
x_perov_left = par_t.dcum0(3);     % dcum is the device thickness

%% Find change in voltage drop across the debye length
%Get index at which the perovskite starts
Debye_left_V_index = find(x==sol_CV_with_ions.par.dcum0(3));
%Get index till which the debye lengths extend to
target=sol_CV_with_ions.par.dcum0(3)+ N_Debye*debye_length; 
temp = abs(target - x);
closest = x(find(temp == min(abs(target - x))));
Debye_right_V_index=find(x==(closest));
%Get difference between the perovskite
deltaV = dfana.deltaVt(sol_CV_with_ions, Debye_left_V_index, Debye_right_V_index);
%% Replicate charges for all time perdids
NA0 = repmat(dev.NA, length(t), 1);
ND0 = repmat(dev.ND, length(t), 1);
Nani0 = repmat(dev.Nani, length(t), 1);
Ncat0 = repmat(dev.Ncat, length(t), 1);
%% Get total, ionic and electronic charge densities across debye length and integrate using trapz to get the
%space charge density
pl = Debye_left_V_index;
pr = Debye_right_V_index;

rho = dfana.calcrho(sol_CV_with_ions,"whole");
Q = e*trapz(x(pl:pr), rho(:, pl:pr), 2);
electronic_rho = -n + p - NA0 + ND0;
EQ = e*trapz(x(pl:pr), electronic_rho(:, pl:pr), 2);
ionic_rho= par.z_a*a + par.z_c*c - par.z_c*Nani0 - par.z_c*Ncat0;
IQ= e*trapz(x(pl:pr), ionic_rho(:, pl:pr), 2);

%% Get capacitance as total charge (per cm2) divided by dV

C_debye_layers = gradient(Q, deltaV);

C_debye_electronic = gradient(EQ, deltaV);

C_debye_ionic = gradient(IQ, deltaV); 


%subplot(2,1,1);
figure(7444)
plot(deltaV, abs(C_debye_layers), deltaV, abs(C_debye_electronic), '-.', deltaV, abs(C_debye_ionic), '--')
legend('Total Capacitance','Electronic Capacitance','Ionic Capacitance')
xlabel('V across debye layers')
ylabel('Capacitances (F/cm^2)')

%% Get capacitance as total charge (per cm2) divided by Vapplied

C_debye_layers = gradient(Q, Vappt);

C_debye_electronic = gradient(EQ, Vappt);

C_debye_ionic = gradient(IQ, Vappt); 


figure(7445)
plot(Vappt, abs(C_debye_layers), Vappt, abs(C_debye_electronic), '-.', Vappt, abs(C_debye_ionic), '--')
legend('Total Capacitance','Electronic Capacitance','Ionic Capacitance')
xlabel('Voltage applied')
ylabel('Capacitances (F/cm^2)')


end
