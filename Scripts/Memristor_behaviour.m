%% Simulate Ag and Iodine reaction by varying vacancy concentration at the interface and introduce
% a surface recombination velocity
initialise_df

%% Define memristor
par_memristor = pc('Input_files/memristor.csv');

%% Get Equilbrium solutions
soleq_memristor = equilibrate(par_memristor);

%% Soleq plot
dfplot.ELnpx(soleq_memristor.ion);
%% Cyclic Voltammogram scan
k_scan = 0.1;
tpoints=200;

Vmax = 6;
Vmin = -6;

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_ion = doCV(soleq_memristor.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
sol_CV_el = doCV(soleq_memristor.el, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
%%No change with sol_JV achieved. Need to find how to vary the surface
%%recombination formula
%% Plot
dfplot.JtotVapp(sol_CV_ion, 0);
hold on
dfplot.JtotVapp(sol_CV_el, 0);
hold off
%set(gca,'YScale','log');
%%
dc=par_memristor.epp0*par_memristor.epp*par_memristor.e;
A=1;
d=par_memristor.d(1);
C=(A*dc)/d;

%%
[J, j, x] = dfana.calcJ(soleq_memristor.ion, "sub");
jn_l = j.n(1);%change to j.n and so on for diffent species (small j means flux)(j solution structure is t_points x x_points)
jn_r = j.n(end);
