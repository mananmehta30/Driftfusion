% Single-layer MAPbICl device, variable workfunctions and ion BC
% 24/03/2021
%
%% Load parameters and customize if necessary
% par = pc('Input_files/1_layer_test.csv');
par = pc('1_layer_MAPI_ITO_Ag.csv');

sc_array = [0,1e-12, 1e-10, 1e-8, 1e-6, 1e-4];
sc_size = length(sc_array);
%% BC range loop

par_temp = par;
% par_temp.Phi_right = -4.6;
% par_temp.Phi_left = -5.0;
for i = 1:sc_size
    par_temp.sc_r = sc_array(i);
    par_temp.sc_l = par_temp.sc_r;
    par_temp = refresh_device(par_temp);
    
    soleq(i) = equilibrate(par_temp);
end

%% Calculate electron only solution
%el_CV = doCV(soleq(1).el, 0, 0, 1.2, -1.2, 1e-1, 2, 241);


%% Plot different BC, both sides, medium scan rate (0.1 Vs-1), two cycles
k = 1e-1;
cycles = 1;

figure()
for i = 1:sc_size
    sol_CV(i) = doCV(soleq(i).ion, 0, 0, 1.2, -1.2, k, cycles, 241);
    dfplot.JtotVapp(sol_CV(i),0)
    hold on
end

dfplot.JtotVapp(el_CV,0)
hold off
% set(gca,'yscale','log')
legentries = cellstr(num2str(sc_array', 'sc=%g'));
legentries{end+1} = 'el only';
legend(legentries)
title(sprintf('%i cycle, scan rate = %g Vs-1, sc on both sides',[cycles,k]))

%% Plot electronic and ionic densities across device, sc both sides
title_arr = strcat('sc = ',string(sc_array));
title_arr(end+1) = 'el only';

sol_CV(end+1) = el_CV;


for i = 1:length(sc_array)
    figure(10+i)
    sgtitle(title_arr(i))
    subplot(2,4,1)
    title('0 V')
    dfplot.npx(sol_CV(i),0/k)
    subplot(2,4,2)
    title('1.1 V')
    dfplot.npx(sol_CV(i),1.1/k)
    subplot(2,4,3)
    title('1.2 V')
    dfplot.npx(sol_CV(i),1.2/k)
    subplot(2,4,4)
    title('-0.5 V')
    dfplot.npx(sol_CV(i), 2.9/k)
    subplot(2,4,5)
    dfplot.acx(sol_CV(i),0/k)
    set(gca,'yscale','log')
    subplot(2,4,6)
    dfplot.acx(sol_CV(i),1.1/k)
    set(gca,'yscale','log')
    subplot(2,4,7)
    dfplot.acx(sol_CV(i),1.2/k)
    set(gca,'yscale','log')
    subplot(2,4,8)
    dfplot.acx(sol_CV(i),2.9/k)
    set(gca,'yscale','log')
end


%% ELx
figure()
dfplot.ELx(sol_CV(1))


%% BC on either side
sc = 1e-8;
sc_l_array = [0,sc,0,sc];
sc_r_array = [0,0,sc,sc];

for i = 1:4
    par_temp.sc_l = sc_l_array(i);
    par_temp.sc_r = sc_r_array(i);
    par_temp = refresh_device(par_temp);
    
    soleq(i) = equilibrate(par_temp);
end

figure()
for i = 1:4
    sol_CV(i) = doCV(soleq(i).ion, 0, 0, 1.2, -1.2, 1e-1, 1, 241);
    
    dfplot.JtotVapp(sol_CV(i),0)
    hold on
end
hold off
legend('no sc','left only','right only','both')

%% Plot carrier densities, sc either side
sol_CV(end+1) = el_CV;
sides = ["none" "left" "right" "both"];
title_arr = strcat('sc = 1e-8, sides: ', sides);
title_arr(end+1) = 'el only';

for i = 1: length(title_arr)
    figure(20+i)
    sgtitle(title_arr(i))
    subplot(2,4,1)
    title('0 V')
    dfplot.npx(sol_CV(i),0)
    subplot(2,4,2)
    title('1.1 V')
    dfplot.npx(sol_CV(i),11)
    subplot(2,4,3)
    title('1.2 V')
    dfplot.npx(sol_CV(i),12)
    subplot(2,4,4)
    title('-0.5 V')
    dfplot.npx(sol_CV(i), 29)
    subplot(2,4,5)
    dfplot.acx(sol_CV(i),0)
    set(gca,'yscale','log')
    subplot(2,4,6)
    dfplot.acx(sol_CV(i),10)
    set(gca,'yscale','log')
    subplot(2,4,7)
    dfplot.acx(sol_CV(i),12)
    set(gca,'yscale','log')
    subplot(2,4,8)
    dfplot.acx(sol_CV(i),29)
    set(gca,'yscale','log')
end

%% Variable scan rates (0.001, 0.01, 0.1, 1 Vs-1)
% sc = 1e-8

sol_CV_0p001 = doCV(soleq(4).ion, 0, 0, 1.2, -1.2, 1e-3, 1, 241);
sol_CV_0p01 = doCV(soleq(4).ion, 0, 0, 1.2, -1.2, 1e-2, 1, 241);
sol_CV_0p1 = doCV(soleq(4).ion, 0, 0, 1.2, -1.2, 1e-1, 1, 241);
sol_CV_1 = doCV(soleq(4).ion, 0, 0, 1.2, -1.2, 1e0, 1, 241);

figure()
dfplot.JtotVapp(sol_CV_0p001,0)
hold on
dfplot.JtotVapp(sol_CV_0p01,0)
dfplot.JtotVapp(sol_CV_0p1,0)
dfplot.JtotVapp(sol_CV_1,0)
dfplot.JtotVapp(el_CV,0)

hold off
legend('0.001 Vs-1', '0.01 Vs-1', '0.1 Vs-1', '1 Vs-1', 'el only, 0.1 Vs-1')
% set(gca,'yscale','log')

%% Plot carrier densities, variable scan rates
sol_CV = [sol_CV_0p001,sol_CV_0p01,sol_CV_0p1,sol_CV_1];
title_arr = strcat('scan rate = ',["0.001","0.01","0.1","1"], 'Vs-1');

for i = 1: length(title_arr)
    figure(20+i)
    sgtitle(title_arr(i))
    subplot(2,4,1)
    title('0 V')
    dfplot.npx(sol_CV(i),0)
    subplot(2,4,2)
    title('1.1 V')
    dfplot.npx(sol_CV(i),11)
    subplot(2,4,3)
    title('1.2 V')
    dfplot.npx(sol_CV(i),12)
    subplot(2,4,4)
    title('-0.5 V')
    dfplot.npx(sol_CV(i), 29)
    subplot(2,4,5)
    dfplot.acx(sol_CV(i),0)
    set(gca,'yscale','log')
    subplot(2,4,6)
    dfplot.acx(sol_CV(i),10)
    set(gca,'yscale','log')
    subplot(2,4,7)
    dfplot.acx(sol_CV(i),12)
    set(gca,'yscale','log')
    subplot(2,4,8)
    dfplot.acx(sol_CV(i),29)
    set(gca,'yscale','log')
end

%% Energy level diagrams
% sc = 5e-9, k = 0.1 Vs-1, 1 loop

sc = 1e-8;
par_temp = pc('Input_files/1_layer_test.csv');
soleq_EL_orig = equilibrate(par_temp);

solCV_EL_orig = doCV(soleq_EL_orig.ion, 0, 0, 1.5, -1.5, 1e-1, 1, 601);

par_temp.sc_l = sc;
par_temp.sc_r = sc;
par_temp = refresh_device(par_temp);

soleq_EL = equilibrate(par_temp);

solCV_EL = doCV(soleq_EL.ion, 0, 0, 1.5, -1.5, 1e-1, 1, 601);

%% ELx comparison
% %  [0,41,121,181,241]
% 
% subplot(2,4,1)
% dfplot.ELx(solCV_EL_orig,0)
% subplot(2,4,1)
% dfplot.ELx(solCV_EL_orig,0)
% title('0 V')
% subplot(2,4,2)
% dfplot.ELx(solCV_EL_orig,15)
% title('+1.5 V')
% subplot(2,4,3)
% dfplot.ELx(solCV_EL_orig,30)
% title('+0 V')
% subplot(2,4,4)
% dfplot.ELx(solCV_EL_orig,45)
% title('-1.5 V')
% subplot(2,4,5)
% dfplot.ELx(solCV_EL,0)
% title('0 V')
% subplot(2,4,6)
% dfplot.ELx(solCV_EL,15)
% title('1.5 V')
% subplot(2,4,7)
% dfplot.ELx(solCV_EL,30)
% title('+0 V')
% subplot(2,4,8)
% dfplot.ELx(solCV_EL,45)
% title('-1.5 V')
% 
% 
% dfplot.JtotVapp(solCV_EL_orig,0)
% hold on
% dfplot.JtotVapp(solCV_EL,0)
% hold off
% legend('no sc',num2str(sc,'sc = %g cm s-1'))
%% change no of ionic species
% 
% par_temp.N_ionic_species = 1;
% par_temp = refresh_device(par_temp);
% 
% sol_eq_1ion = equilibrate(par_temp);
% 
% sol_CV_1ion = doCV(sol_eq_1ion.ion,0, 0, 1.2, -1.2, 1e-1, 1, 241);

%%
% dfplot.JtotVapp(sol_CV_1ion,0)
% hold on
% dfplot.JtotVapp(sol_CV_0p1,0)
% hold off
% legend('1 ionic species','2 ionic species')
