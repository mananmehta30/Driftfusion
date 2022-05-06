function [DeltaV, Q, C] = capacitance_ana_PC(sol, pervsk_layer_no)

[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

%% Calculate Debye length from Schottky defect density and perovskite high frequency dielectric constant
% I have decided in this version to use the midpoint of the perovskite for
% the charge integration step
e = par.e;
V_T = par.kB*par.T;                     % Thermal votlage
epp_pvsk = e*par.epp0*par.epp(3);       % Perovskite absolute dielectric constant
N0 = par.Ncat(3);
L_D = sqrt((epp_pvsk*V_T)/(e*N0));      % Debye length

%% Get point indices for various locations
pcum0 = par.pcum0;
dcum0 = par.dcum0;

pl_SCR1 = pcum0(pervsk_layer_no + 1);
pr_SCR1 = round((pcum0(pervsk_layer_no + 1) + pcum0(pervsk_layer_no + 2))/2);   % use mid point of pvsk as first approx - would be better to use point at which sign of charge changes
pl_SCR2 = round((pcum0(pervsk_layer_no + 1) + pcum0(pervsk_layer_no + 2))/2);
pr_SCR2 = pcum0(pervsk_layer_no + 2);

%% Calculate volumetric charge components (cm-3)
NA = repmat(dev.NA, length(t), 1);
ND = repmat(dev.ND, length(t), 1);
Nani = repmat(dev.Nani, length(t), 1);
Ncat = repmat(dev.Ncat, length(t), 1);

rho = dfana.calcrho(sol, "whole");
rho_el = -n + p - NA + ND;
rho_ion = par.z_a*a + par.z_c*c - par.z_c*Nani - par.z_c*Ncat;

%% Get areal charge densities (C cm-2) for different SCRs
Q.SCR1 = e*trapz(x(pl_SCR1:pr_SCR1), rho(:, pl_SCR1:pr_SCR1), 2);
Q.el_SCR1 = e*trapz(x(pl_SCR1:pr_SCR1), rho_el(:, pl_SCR1:pr_SCR1), 2);
Q.ion_SCR1 = e*trapz(x(pl_SCR1:pr_SCR1), rho_ion(:, pl_SCR1:pr_SCR1), 2);

Q.SCR2 = e*trapz(x(pl_SCR2:pr_SCR2), rho(:, pl_SCR2:pr_SCR2), 2);
Q.el_SCR2 = e*trapz(x(pl_SCR2:pr_SCR2), rho_el(:, pl_SCR2:pr_SCR2), 2);
Q.ion_SCR2 = e*trapz(x(pl_SCR2:pr_SCR2), rho_ion(:, pl_SCR2:pr_SCR2), 2);

Q.pvsk = e*trapz(x(pl_SCR1:pr_SCR2), abs(rho(:, pl_SCR1:pr_SCR2)), 2);

%% Get Voltage changes across the SCRs
DeltaV.SCR1 = V(:, pr_SCR1) - V(:, pl_SCR1);
DeltaV.SCR2 = V(:, pr_SCR2) - V(:, pl_SCR2);
DeltaV.pvsk = V(:, pr_SCR2) - V(:, pl_SCR1);

%% Calculate differential capacitance dQ/dV
C.SCR1 = abs(gradient(Q.SCR1, DeltaV.SCR1));
C.el_SCR1 = abs(gradient(Q.el_SCR1, DeltaV.SCR1));
C.ion_SCR1 = abs(gradient(Q.ion_SCR1, DeltaV.SCR1));

C.SCR2 = abs(gradient(Q.SCR2, DeltaV.SCR2));
C.el_SCR2 = abs(gradient(Q.el_SCR2, DeltaV.SCR2));
C.ion_SCR2 = abs(gradient(Q.ion_SCR2, DeltaV.SCR2));

%% Total capacitance of perovskite as series sum
C.tot = (1./C.SCR1 + 1./C.SCR2).^-1;
C.el = (1./C.el_SCR1 + 1./C.el_SCR2).^-1;
C.ion = (1./C.ion_SCR1 + 1./C.ion_SCR2).^-1;

%% Get total applied voltage
Vapp = dfana.calcVapp(sol);

%% Plot the voltages
figure(597)
plot(t, Vapp, t, DeltaV.pvsk, t, DeltaV.SCR1, t, DeltaV.SCR2)
xlabel("Time (s)")
ylabel("Voltage (V)")
legend("V_{app}", "\Delta V_{pvsk}", "\Delta V_{SCR1}", "\Delta V_{SCR2}")

%% Plot the charges
figure(598)
plot(t, Q.pvsk, t, Q.SCR1, t, Q.SCR2)
xlabel("Time (s)")
ylabel("Charge density (Ccm^{-2})")
legend("|Q_{SCR1}| + |Q_{SCR2}|", "Q_{SCR1}", "Q_{SCR2}")

%% Plot the capacitances as function of Vapp
figure(599)
semilogy(Vapp, C.tot, Vapp, C.SCR1, '--', Vapp, C.SCR2, '--')
ylabel("Capacitance (Fcm^{-2})")
xlabel("Gate voltage (V)")
legend("Perovskite total", "SCR1", "SCR2")

figure(600)
semilogy(Vapp, C.tot, Vapp, C.ion, '-.', Vapp, C.ion_SCR1, '--', Vapp, C.ion_SCR2, '--')
ylabel("Capacitance (Fcm^{-2})")
xlabel("Gate voltage (V)")
legend("Perovskite total", "Ionic", "SCR1 - Ionic", "SCR2 - Ionic")

figure(601)
semilogy(Vapp, C.tot, Vapp, C.el, '-.', Vapp, C.el_SCR1, '--', Vapp, C.el_SCR2, '--')
ylabel("Capacitance (Fcm^{-2})")
xlabel("Gate voltage (V)")
legend("Perovskite total", "Electronic", "SCR1", "SCR1 - Electronic", "SCR2", "SCR2 - Electronic")    

end