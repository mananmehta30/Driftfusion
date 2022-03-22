function [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV)

par = sol_CV.par;

%% Calculate conductivity
[sigma_n, sigma_p] = dfana.calc_conductivity(sol_CV);

%% Debye length Calculation
e = par.e;
V_T = par.kB*par.T;                     % Thermal votlage
epp_pvsk = e*par.epp0*par.epp(3);       % Perovskite absolute dielectric constant
N0 = par.Ncat(3);
L_D = sqrt((epp_pvsk*V_T)/(e*N0));      % Deby width [cm]
N_Debye = 5;                            % Number of Debye lengths to average electron density over
x_perov_left = sol_CV.par.dcum0(3);     % dcum is the device thickness
x_perov_right = sol_CV.par.dcum0(4);
x = sol_CV.x;
t = sol_CV.t;
Vappt = dfana.calcVapp(sol_CV);
%% Find mean conductivity
sigma_n_bar = sigma_n(:, sol_CV.par.pcum0(3) +1); %
sigma_p_bar = sigma_p(:, sol_CV.par.pcum0(3) +1);

%% Find peak conductivity for applied bias
pp_Vmax = find(Vappt == max(Vappt));      %% pp = point position
pp_Vmin = find(Vappt == min(Vappt));      %% pp = point position

%% Max conductivity
sigma_n_bar_Vpeak = sigma_n_bar(pp_Vmax);
sigma_p_bar_Vpeak = sigma_p_bar(pp_Vmax);

end