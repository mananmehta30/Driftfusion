%% Code pupose
% Code is supposed to take different values of doping in different columns and different left 
% electrode workfunction values in the rows and give max conductivity achieved
%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');

par = par_alox;     % Create temporary parameters object for overwriting parameters in loop

%% Initialise the parameter arrays
Ncat_array = logspace(16, 19, 10);

thickness_array = 0.000005:0.000005:0.00003;


%% while
for i = 1:length(Ncat_array)
    
    par.Ncat(:) = Ncat_array(i);
    par.Nani(:) = Ncat_array(i);
    
    disp(['Cation density = ', num2str(Ncat_array(i)), ' cm^-3']);
    for j = 1:length(thickness_array) %loop to run for insulator thickness
        

        par.d(1) = thickness_array(j);

        disp(['Insulator Thickness = ', num2str(thickness_array(j)), 'cm']);

 

        
        par = refresh_device(par);      % This line is required to rebuild various arrays used DF
        
        %% Find equilibrium
        soleq(i, j) = equilibrate(par);
        dfplot.acx(soleq(i, j).ion)
        
        %% Current-voltage scan
        k_scan = 0.001;
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
sigma_n_barM = zeros(length(Ncat_array), length(thickness_array), length(sol_CV(1,1).t)); 
sigma_p_barM = zeros(length(Ncat_array), length(thickness_array), length(sol_CV(1,1).t)); 
sigma_n_bar_VpeakM = zeros(length(Ncat_array), length(thickness_array)); 
sigma_p_bar_VpeakM = zeros(length(Ncat_array), length(thickness_array)); 

for i = 1:length(Ncat_array)
    for j = 1:length(thickness_array)
        [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV(i,j));
        sigma_n_barM(i,j,:) = sigma_n_bar;
        sigma_p_barM(i,j,:) = sigma_p_bar;
        sigma_n_bar_VpeakM(i,j) = sigma_n_bar_Vpeak;
        sigma_p_bar_VpeakM(i,j) = sigma_p_bar_Vpeak;
    end
end

%% Plots

%% Conductivity vs Insulator Thickness
for i = 1:length(Ncat_array)
    figure(100)
    semilogy(thickness_array, sigma_n_bar_VpeakM(i, :))
    hold on
    xlabel('Insulator Thickness [cm]')
    ylabel('Peak electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
end  

for i = 1:length(Ncat_array)
    figure(101)
    semilogy(thickness_array, sigma_p_bar_VpeakM(i, :))
    hold on
    xlabel('Insulator Thickness [cm]]')
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


insulator_thickness_index = 2;
for i=1:241
    for j=1:length(Ncat_array)
    conductivity(j,i)=sigma_n_barM(j, insulator_thickness_index,i);

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
for j = 1:length(thickness_array)
    figure(201)
    semilogy(Vappt, squeeze(sigma_n_barM(3, j, :)))
    legstr_n2{j} = ['Thickness = ', num2str(thickness_array(j)),' ','cm'];
    hold on
end

for j = 1:length(thickness_array)
    figure(202)
    semilogy(Vappt, squeeze(sigma_p_barM(3, j, :)))
    legstr_p2{j} = ['Thickness =', num2str(thickness_array(j)),' ','cm'];
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


insultor_thickness_index = 3;

legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)

    n_int = sol_CV(i, insultor_thickness_index).u(:, par.pcum0(3), 2);

    figure(203)
    semilogy(Vappt, n_int)
    legstr_n3{i} = ['Ncat =', num2str(Z(i)),' cm-3'];
    hold on
end

for i = 1:length(Ncat_array)

    p_int = sol_CV(i, insultor_thickness_index).u(:, par.pcum0(3), 3);

    figure(204)
    semilogy(Vappt, p_int)
    legstr_p3{i} = ['Ncat =', num2str(Z(i)),' cm-3'];
    hold on
end

figure(203)
xlabel('Voltage [V]')
ylabel('Electron density at the interface (cm-3)')
legend(legstr_n3)
hold off

figure(204)
xlabel('Voltage [V]')
ylabel('Hole density interface (cm-3)')
legend(legstr_p3)
hold off

%% Plot electron and hole profiles at Vmax as a function of position

legstr_npx = {'', '', ''};
for i = 1:length(Ncat_array)

    dfplot.npx(sol_CV(i, insultor_thickness_index), Vmax/k_scan);% Vmax/k_scan)

    dfplot.npx(sol_CV(i, insulator_thickness_index), Vmax/k_scan);% Vmax/k_scan)

    legstr_npx{2*i-1 + 3} = ['n, Ncat =', num2str(Z(i)),' cm-3'];
    legstr_npx{2*i + 3} = ['p, Ncat =', num2str(Z(i)),' cm-3'];
    hold on
end
legend(legstr_npx)
ylim([1e-1, 1e12])

%% Plot cation and anion densities at Vmax as a function of position




%dfplot.acx(sol_CV(3, 9), 3*(Vmax/k_scan)
%% Plot potential as a function position

insulator_thickness_index = 1;
legstr_Vx = {'dielectric', 'interface', 'perovskite'};
for i = 1:length(Ncat_array)
    dfplot.Vx(sol_CV(i, insulator_thickness_index), Vmax/k_scan);%Vmax/k_scan)
     Z(i)=round(Ncat_array(i),3,'significant');
    legstr_Vx{i + 3} = ['Ncat =', num2str(Z(i)),' cm-3'];
    hold on
end
legend(legstr_Vx)
%ylim([1e-1, 1e12])



%% Plot average conductivity
% figure(201)
% plot(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
% xlabel('Voltage [V]')
% ylabel('Average channel conductivity [Linear]')
% legend('Electron', 'Hole')

%% Plot average conductivity
% figure(202)
% semilogy(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
% xlabel('Voltage [V]')
% ylabel('Average bulk conductivity [Semilog]')
% legend('Electron', 'Hole')

%% Plot average conductivity
% figure(203)
% plot(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
% xlabel('Voltage [V]')
% ylabel('Average bulk conductivity [Linear]')
% legend('Electron', 'Hole')

%% Plot average conductivity
% figure(204)
% semilogy(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
% xlabel('Voltage [V]')
% ylabel('Average entire conductivity [Semilog]')
% legend('Electron', 'Hole')

%% Plot average conductivity
% figure(205)
% plot(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
% xlabel('Voltage [V]')
% ylabel('Average entire conductivity [Linear]')
% legend('Electron', 'Hole')


%% Plot ionic concentration at interface as function Vapp for different ionic densities


%Plot similar for ions (instead of n_int put cat_int)

insultor_thickness_index = 4;


legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)

    cat_int = sol_CV(i, insultor_thickness_index).u(:, par.pcum0(3)+1,4);
    logcat_int=log10(cat_int);
    figure(703)
    semilogy(Vappt,cat_int)
    Z(i)=round(Ncat_array(i),3,'significant');
    legstr_n3{i} = ['Ncat =', num2str(Z(i))];
    hold on
end
figure(703)
xlabel('Voltage [V]')
ylabel('Cation Concentration (cm-3)')
legend(legstr_n3)
hold off
%% Modulability Ions


Ncat_index = 4;
legstr_n3 =[];
legstr_p3 =[];
for i = 1:length(thickness_array)
    cat_int = sol_CV(Ncat_index, i).u(:, par.pcum0(3)+1,4);

    cat_int_log=log10(cat_int);%log(c)

           cation_modulability=gradient(cat_int_log,Vappt);%dlog(n)/dV
            
            cation_modulability_factor(i)= abs(cation_modulability(1));
   
end

figure(1112)
scatter(thickness_array, cation_modulability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','linear')

%xlim([1e15 1e20])
%ylim([4.6 8.5])

xlabel('Insulator Thickness')
ylabel('Cation Modulability Factor (m_V_g)')
box on

%% Electon concentration Modulability vs Cation Concentration
insultor_thickness_index=2;
for i = 1:length(Ncat_array)
    
        built_in_potential=0;
          n_int = sol_CV(i, insultor_thickness_index).u(:, par.pcum0(3), 2);
          log_n=log10(n_int);%log(n)
           n_modulability=gradient(log_n,Vappt);%dlog(n)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            n_modulability_factor(i)= n_modulability(I);
    
end

figure(1112)
scatter(Ncat_array, n_modulability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
%ylim([4.6 8.5])

xlabel('Cation concentration')
ylabel('Electron Modulability Factor (m_V_g)')
legend('Modulatability Factor')
box on







%% Electon concentration Modulability vs Insulator Thickness
Ncat_index=10;
for i = 1:length(thickness_array)
    
        built_in_potential=0;
          nnn_int = sol_CV(Ncat_index, i).u(:, par.pcum0(3), 2);
          log_nnn=log10(nnn_int);%log(n)
           nnn_modulability=gradient(log_nnn,Vappt);%dlog(n)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            nnn_modulability_factor(i)= nnn_modulability(1);
    
end

figure(1155)
scatter(thickness_array, nnn_modulability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','linear')

xlim([0 3.5e-5])
ylim([0 0.55])

xlabel('Insulator Thickess [cm]')
ylabel('Electron Modulability Factor [m_V_g] ')
legend('Modulatability Factor for Cation Density = 1e19')

box on









%% Electon conductivity Modulability vs Cation Concentration
insultor_thickness_index=1;
for i = 1:length(Ncat_array)
    
        
          sigma_n_int = sigma_n_bar;
          log_sigma_n=log10(sigma_n_int);%log(n)
           sigma_n_modulability=gradient(log_sigma_n,Vappt);%dlog(sigma)/dV
            sigma_n_modulability_factor(i)= sigma_n_modulability(1);
    
end

figure(2222)
scatter(Ncat_array, sigma_n_modulability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
%ylim([4.6 8.5])

xlabel('Cation concentration')
ylabel('Electron Conductivity Modulability Factor (m_V_g)')
box on
%% Hole concentration Modulability vs Cation Concentration
insultor_thickness_index=1;
for i = 1:length(Ncat_array)

          p_int = sol_CV(i, insultor_thickness_index).u(:, par.pcum0(3), 3);
          log_p=log10(p_int);%log(n)
           p_modulability=gradient(log_p,Vappt);%dlog(p)/dV
            p_modulability_factor(i)= p_modulability(I);
    
end

figure(1113)
scatter(Ncat_array, p_modulability_factor,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
%ylim([4.6 8.5])

xlabel('Cation concentration')
ylabel('Hole Modulability Factor (m_V_g)')
box on
%% Electron Modulability Contour



for i = 1:length(Ncat_array)
    for j=1:length(thickness_array)
        built_in_potential=0;
          n_int = sol_CV(i, j).u(:, par.pcum0(3), 2);
          log_n=log10(n_int);%log(n)
            n_modulability=gradient(log_n,Vappt);%dlog(n)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            n_modulability_factor_contour(i,j)= n_modulability(I);
    end 
end

x=thickness_array;
y=Ncat_array;
z=n_modulability_factor_contour;
z_log=log10(z);
figure(1)
surf(x,y,z);
set(gca,'ZScale','linear')
xlabel('Insulator thickness'), ylabel('Cation Concentration'), zlabel('Modulability factor')
set(gca,'YScale','log')
box on
%% Hole Modulability Contour



for i = 1:length(Ncat_array)
    for j=1:length(thickness_array)
        built_in_potential=0;
          p_int = sol_CV(i, j).u(:, par.pcum0(3), 3);
          log_p=log10(p_int);%log(n)
            p_modulability=gradient(log_p,Vappt);%dlog(p)/dV
            target=built_in_potential; 
             temp=abs(target-Vappt);
             [M,I] = min(temp);
            p_modulability_factor_contour(i,j)= p_modulability(I);
    end 
end
figure(2)
x2=thickness_array;
y2=Ncat_array;
z2=p_modulability_factor_contour;
z2_log=log10(z2);
surf(x2,y2,z2);

set(gca,'ZScale','linear')
xlabel('Workfunction'), ylabel('Cation Concentration'), zlabel('Modulability factor')
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
