% Single carrier device test
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code
%% Load parameters
par = pc('Input_files/1_layer_single_carrier_Ag.csv');

%% Run to equilibrium
soleq_single_carrier = equilibrate(par);

%% Initial ennergy level plot
dfplot.ELx(soleq_single_carrier.ion);
% hold on
% scatter([0,par.d*1e7],[par.Phi_left,par.Phi_right])
% hold off
%% Do Cyclic Voltammograms (CV) at 1, 10, and 100 Vs-1 for 4 cycles each
% doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_2cyc = doCV(soleq_single_carrier.ion, 0, 0, 0.6, -0.2, 1e-2, 2, 161);
sol_CV_4cyc = doCV(soleq_single_carrier.ion, 0, 0, 0.6, -0.2, 1e-2, 4, 321);
sol_CV_6cyc = doCV(soleq_single_carrier.ion, 0, 0, 0.6, -0.2, 1e-2, 6, 481);

%% Plot CV using currents from left-hand boundary (x=0)
dfplot.JtotVapp(sol_CV_2cyc, 0);
hold on
dfplot.JtotVapp(sol_CV_4cyc, 0);
hold on
dfplot.JtotVapp(sol_CV_6cyc, 0);
legend('2 cycles', '4 cycles','6 cycles')
hold off
%% Parameters exploration loop
% sc_l_array = [0, 1e-12, 1e-10, 1e-8, 1e-6];
% par_temp = par;
% 
% for i=1:length(sc_l_array)
%     par_temp.sc_l = sc_l_array(i);
%     par_temp = refresh_device(par_temp);
%     
%     soleq(i) = equilibrate(par_temp);
%     sol_CV(i) = doCV(soleq(i).ion, 0, 0, 0.6, -0.6, 1e-2, 2, 241);
% end

%% Plot 
% for i=1:length(sc_l_array)
%     dfplot.JtotVapp(sol_CV(i), 0);
%     hold on
% end
% hold off
% legend(cellstr(num2str(sc_l_array', 'sc_l=%-d')));

%% Make movie
% makemovie(sol_CV_k0p01, @dfplot.acx, 0, 0, 'ELx', 1, 0)