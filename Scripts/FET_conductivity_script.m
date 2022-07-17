%% Code pupose
%% These results can prove the scan rate or time scale thingy

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');
par = par_alox;     % Create temporary parameters object for overwriting parameters in loop

%% Initialise the parameter arrays
%Ncat_array = logspace(16, 19, 4);
%Ncat_array=[1e16,5e16,1e17,5e17,1e18,5e18,1e19];
mobility_arrays = [0,1e-10,1e-9,1e-7,1e-5];
Ncat_array=logspace(16,19,10);
par.Phi_left=-4.9;
par.Phi_right=-4.9;
%% No of ionic species
par.N_ionic_species=1;
%% while
for i = 1:length(Ncat_array)
    
    par.Ncat(:) = Ncat_array(i);
    par.Nani(:) = Ncat_array(i);
    
    disp(['Cation density = ', num2str(Ncat_array(i)), ' cm^-3']);
    for j = 1:length(mobility_arrays) %loop to run for different electrode workfunction
        
        par.mu_c(3) = mobility_arrays(j);
        disp(['Mobility of cation = ', num2str(mobility_arrays(j)), ' cm2/Vs']);
        
        par = refresh_device(par);      % This line is required to rebuild various arrays used DF
        
        %% Find equilibrium
        soleq(i, j) = equilibrate(par);
        dfplot.acx(soleq(i, j).ion)
        
        %% Current-voltage scan
        k_scan = 0.1;
        Vmax = 1.2;
        Vmin = -1.2;
        
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        sol_CV(i, j) = doCV(soleq(i, j).ion, 0, 0, Vmax, Vmin, k_scan, 1, 241);
        %% Plot Vapp vs time
        % dfplot.Vappt(sol_CV)
        
        %% Plot JV scan
        %dfplot.JtotVapp(sol_CV, 0);
        %set(gca,'YScale','log')
        
        %% Plot anion and cation densities
        %dfplot.acx(sol_CV, 1/k_scan*[0:Vmax/3:Vmax]);
        
        %% Plot electron and hole profiles
        %dfplot.npx(sol_CV, 1/k_scan*[0:Vmax/3:Vmax]);
    end   
end

%% Analysis
Vappt = dfana.calcVapp(sol_CV(1,1));
% Preallocation
sigma_n_barM = zeros(length(Ncat_array), length(mobility_arrays), length(sol_CV(1,1).t)); 
sigma_p_barM = zeros(length(Ncat_array), length(mobility_arrays), length(sol_CV(1,1).t)); 
sigma_n_bar_VpeakM = zeros(length(Ncat_array), length(mobility_arrays)); 
sigma_p_bar_VpeakM = zeros(length(Ncat_array), length(mobility_arrays)); 

for i = 1:length(Ncat_array)
    for j = 1:length(mobility_arrays)
        [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV(i,j));
        sigma_n_barM(i,j,:) = sigma_n_bar;
        sigma_p_barM(i,j,:) = sigma_p_bar;
        sigma_n_bar_VpeakM(i,j) = sigma_n_bar_Vpeak;
        n_bar_VpeakM(i,j)= (sigma_n_bar_VpeakM(i,j)./par.e)./par.mu_n(3);
        sigma_p_bar_VpeakM(i,j) = sigma_p_bar_Vpeak;
        p_bar_VpeakM(i,j)= (sigma_p_bar_VpeakM(i,j)./par.e)./par.mu_p(3);
    end
end

%% Plots

%% Conductivity vs Work Function
for i = 1:length(Ncat_array)
    figure(100)
    semilogy(mobility_arrays, sigma_n_bar_VpeakM(i, :))
    hold on
    xlabel('Mobility')
    ylabel('Peak electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
end  

for i = 1:length(Ncat_array)
    figure(101)
    semilogy(mobility_arrays, sigma_p_bar_VpeakM(i, :))
    hold on
    xlabel('Mobility')
    ylabel('Peak hole conductivity [S cm-1]')
    legstr_p{i} = ['Ncat =', num2str(Ncat_array(i))];
end  
figure(100)
legend(legstr_n)
hold off
figure(101)
legend(legstr_p)
hold off


%% Conductivity vs Applied Voltage for different ionic densities

mobility_index = 1;
for i=1:241
    for j=1:length(Ncat_array)
    conductivity(j,i)=sigma_n_barM(j, mobility_index,i);
    end
end

for i = 1:length(Ncat_array)
    figure(1000)
    semilogy(Vappt, conductivity(i,:));
    hold on
    xlabel('Applied Voltage [V]')
    ylabel('Electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
end  
figure(1000)
legend(legstr_n)
hold off

%% Plot average conductivity
for j = 1:length(mobility_arrays)
    figure(201)
    semilogy(Vappt, squeeze(sigma_n_barM(3, j, :)))
    legstr_n2{j} = ['Mobility =', num2str(mobility_arrays(j))];
    hold on
end

for j = 1:length(mobility_arrays)
    figure(202)
    semilogy(Vappt, squeeze(sigma_p_barM(3, j, :)))
    legstr_p2{j} = ['Mobility=', num2str(mobility_arrays(j))];
    hold on
end
%%check how to write siemens properly
figure(201)
xlabel('Voltage [V]')
ylabel('Average electron conductivity [S cm-1]')
legend(legstr_n2)
hold off

figure(202)
xlabel('Voltage [V]')
ylabel('Average hole conductivity [S cm-1]')
legend(legstr_p2)
hold off


%% Plot carrier concentration at interface as function Vapp for different ion densities

mobility_index = 5;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    n_int = sol_CV(i, mobility_index).u(:, par.pcum0(3), 2);
    figure(203)
    semilogy(Vappt, n_int)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

for i = 1:length(Ncat_array)
    p_int = sol_CV(i, mobility_index).u(:, par.pcum0(3), 3);
    figure(204)
    semilogy(Vappt, p_int)
    legstr_p3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

figure(203)
xlabel('Applied Voltage [V]')
ylabel('Electron density at the interface (cm-3)')
legend(legstr_n3)
title(legend,'Cation defect density (cm-3)')
hold off

figure(204)
xlabel('Applied Voltage [V]')
ylabel('Hole density interface (cm-3)')
legend(legstr_p3)
title(legend,'Cation defect density (cm-3)')
hold off

%% Plot electron and hole profiles at Vmax as a function of position
mobility_index=7;
legstr_npx = {'', '', ''};
for i = 1:length(Ncat_array)
    dfplot.npx(sol_CV(i, mobility_index), Vmax/k_scan);% Vmax/k_scan)
    legstr_npx{2*i-1 + 3} = ['n, Ncat =', num2str(Ncat_array(i))];
    legstr_npx{2*i + 3} = ['p, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_npx)
ylim([1e-1, 1e12])


%% Plot potential as a function position
mobility_index = 7;
legstr_Vx = {'dielectric', 'interface', 'perovskite'};
for i = 1:length(Ncat_array)
    dfplot.Vx(sol_CV(i, mobility_index), Vmax/k_scan);%Vmax/k_scan)
    legstr_Vx{i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_Vx)
%ylim([1e-1, 1e12])




%% Electon concentration Modulatability vs Cation Concentration
mobility_index=7;
for i = 1:length(Ncat_array)
    
       
          n_int = sol_CV(i, mobility_index).u(:, par.pcum0(3), 2);
          log_n=log10(n_int);%log(n)
           n_Modulatability=gradient(log_n,Vappt);%dlog(n)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            n_Modulatability_factor(i)= n_Modulatability(1);
    
end

figure(1112)
scatter(Ncat_array, n_Modulatability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
%ylim([4.6 8.5])
legend('Modulatability factor (m_V_g)')
xlabel('Ionic concentration (cm-3)')
ylabel('Electron Modulatability Factor (m_V_g)')
box on

%% Calculate manually
mobility_index=7;
for i = 1:length(Ncat_array)
    
        built_in_potential=par.Phi_right-mobility_arrays(mobility_index);
          n_int = sol_CV(i, mobility_index).u(:, par.pcum0(3), 2);
          sigma_nn_int= par.e.*par.mu_n(3).*n_int;

         
          log_nn=log10(sigma_nn_int);%log(n)
           nn_Modulatability=gradient(log_nn,Vappt);%dlog(n)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            nn_Modulatability_factor(i)= nn_Modulatability(I);
    
end
figure(1112)
scatter(Ncat_array, nn_Modulatability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
%ylim([4.6 8.5])

xlabel('Cation concentration (cm-3)')
ylabel('Electron Conductivity Modulatability Factor (m_V_g)')
box on
%% Electon conductivity Modulatability vs Cation Concentration
mobility_index=7;
for i = 1:length(Ncat_array)
    
        built_in_potential=par.Phi_right-mobility_arrays(mobility_index);
          sigma_n_int = sigma_n_bar;
          log_sigma_n=log10(sigma_n_int);%log(n)
           sigma_n_Modulatability=gradient(log_sigma_n,Vappt);%dlog(sigma)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            sigma_n_Modulatability_factor(i)= sigma_n_Modulatability(I);
    
end

figure(2222)
scatter(Ncat_array, sigma_n_Modulatability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
%ylim([4.6 8.5])

xlabel('Cation concentration')
ylabel('Electron Conductivity Modulatability Factor (m_V_g)')
box on
%% Hole concentration Modulatability vs Cation Concentration
mobility_index=7;
for i = 1:length(Ncat_array)
    
        built_in_potential=par.Phi_right-mobility_arrays(mobility_index);
          p_int = sol_CV(i, mobility_index).u(:, par.pcum0(3), 3);
          log_p=log10(p_int);%log(n)
           p_Modulatability=gradient(log_p,Vappt);%dlog(p)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            p_Modulatability_factor(i)= p_Modulatability(I);
    
end

figure(1113)
scatter(Ncat_array, p_Modulatability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
%ylim([4.6 8.5])

xlabel('Cation concentration')
ylabel('Hole Modulatability Factor (m_V_g)')
box on
%% Electron Modulatability Contour



for i = 1:length(Ncat_array)
    for j=1:length(mobility_arrays)
        built_in_potential=par.Phi_right-mobility_arrays(j);
          n_int = sol_CV(i, j).u(:, par.pcum0(3), 2);
          log_n=log10(n_int);%log(n)
            n_Modulatability=gradient(log_n,Vappt);%dlog(n)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            n_Modulatability_factor_contour(i,j)= n_Modulatability(I);
    end 
end

x=mobility_arrays;
y=Ncat_array;
z=n_Modulatability_factor_contour;
z_log=log10(z);
figure(1)
surf(x,y,z);
set(gca,'ZScale','linear')
xlabel('Workfunction'), ylabel('Cation Concentration'), zlabel('Modulatability factor')
set(gca,'YScale','log')
box on
figure(444)
contour(x,y,z)
set(gca,'ZScale','linear')
xlabel('Workfunction'), ylabel('Cation Concentration'), zlabel('Modulatability factor')
set(gca,'YScale','log')
box on

%% Hole Modulatability Contour



for i = 1:length(Ncat_array)
    for j=1:length(mobility_arrays)
        built_in_potential=par.Phi_right-mobility_arrays(j);
          p_int = sol_CV(i, j).u(:, par.pcum0(3), 3);
          log_p=log10(p_int);%log(n)
            p_Modulatability=gradient(log_p,Vappt);%dlog(p)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            p_Modulatability_factor_contour(i,j)= p_Modulatability(I);
    end 
end
figure(2)
x2=mobility_arrays;
y2=Ncat_array;
z2=p_Modulatability_factor_contour;
z2_log=log10(z2);
surf(x2,y2,z2);

set(gca,'ZScale','linear')
xlabel('Workfunction'), ylabel('Cation Concentration'), zlabel('Modulatability factor')
set(gca,'YScale','log')
box on
%% Conductivity profiles
% So systematically you could look at the following.

% 1)Electrode workfunctions
% 2)Ion density 
% 3)1 ion, opposite charge (i.e. mobile anions)
% 4)2 ions
% 5)Different ion densities

% Ideally you would set some of these parameter explorations up as loops and extract peak conductivity then plot 
% on a contour plot with x = Ion density, y = Electrode workfunctions, z = peak conductivity for example.

%% Modulatability Ions

% v_built_in=par.Phi_left-par.Phi_right;
% mobility_index = 3;
% legstr_n3 =[];
% legstr_p3 =[];
% for i = 1:length(Ncat_array)
%     cat_int = sol_CV(i, mobility_index).u(:, par.pcum0(3)+1,4);
%     cat_int_log=log10(cat_int);
%     figure(112)
%     plot(Vappt, cat_int_log)
%     legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
%     hold on
% end
% figure(112)
% xlabel('Voltage [V]')
% ylabel('Cation Concentration (cm-3)')
% 
% hold off

%% Plot average conductivity
% figure(200)
% semilogy(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
% xlabel('Voltage [V]')
% ylabel('Average channel conductivity [Semilog]')
% legend('Electron', 'Hole')
% 
% %% Plot average conductivity
% % figure(201)
% % plot(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
% xlabel('Voltage [V]')
% % ylabel('Average channel conductivity [Linear]')
% % legend('Electron', 'Hole')
% 
% %% Plot average conductivity
% figure(202)
% semilogy(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
% xlabel('Voltage [V]')
% ylabel('Average bulk conductivity [Semilog]')
% legend('Electron', 'Hole')
% 
% %% Plot average conductivity
% % figure(203)
% % plot(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
% % xlabel('Voltage [V]')
% % ylabel('Average bulk conductivity [Linear]')
% % legend('Electron', 'Hole')
% 
% %% Plot average conductivity
% figure(204)
% semilogy(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
% xlabel('Voltage [V]')
% ylabel('Average entire conductivity [Semilog]')
% legend('Electron', 'Hole')
% 
% %% Plot average conductivity
% % figure(205)
% % plot(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
% % xlabel('Voltage [V]')
% % ylabel('Average entire conductivity [Linear]')
% % legend('Electron', 'Hole')
% 
% 


%% Plot ionic concentration at interface as function Vapp for different ion densities


%Plot similar for ions (instead of n_int put cat_int)
% mobility_index = 3;
% legstr_n3 =[];
% legstr_p3 =[];
% 
% for i = 1:length(Ncat_array)
%     cat_int = sol_CV(i, mobility_index).u(:, par.pcum0(3)+1,4);
%     logcat_int=log10(cat_int);
%     figure(703)
%     semilogy(Vappt,cat_int)
%     legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
%     hold on
% end
% figure(703)
% xlabel('Voltage [V]')
% ylabel('Cation Concentration (cm-3)')
% legend(legstr_n3)
% hold off

%% Plot cation and anion densities at Vmax as a function of position




%dfplot.acx(sol_CV(3, 9), 3*(Vmax/k_scan)