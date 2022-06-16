%% Code pupose

%%% Make table of different reactions/theories to explain resistive
%%% switching. Send to all
%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');

par = par_alox;     % Create temporary parameters object for overwriting parameters in loop

%% Initialise the parameter arrays
Ncat_array = logspace(16, 19, 4);
workfunction_LHS = -5.5:0.1:-4.2;

%% while
for i = 1:length(Ncat_array)
    
    par.Ncat(:) = Ncat_array(i);
    par.Nani(:) = Ncat_array(i);
    
    disp(['Cation density = ', num2str(Ncat_array(i)), ' cm^-3']);
    for j = 1:length(workfunction_LHS) %loop to run for different electrode workfunction
        
        par.Phi_left = workfunction_LHS(j);
        disp(['LHS electrode workfunction = ', num2str(workfunction_LHS(j)), ' eV']);
        
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
sigma_n_barM = zeros(length(Ncat_array), length(workfunction_LHS), length(sol_CV(1,1).t)); 
sigma_p_barM = zeros(length(Ncat_array), length(workfunction_LHS), length(sol_CV(1,1).t)); 
sigma_n_bar_VpeakM = zeros(length(Ncat_array), length(workfunction_LHS)); 
sigma_p_bar_VpeakM = zeros(length(Ncat_array), length(workfunction_LHS)); 

for i = 1:length(Ncat_array)
    for j = 1:length(workfunction_LHS)
        [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV(i,j));
        sigma_n_barM(i,j,:) = sigma_n_bar;
        sigma_p_barM(i,j,:) = sigma_p_bar;
        sigma_n_bar_VpeakM(i,j) = sigma_n_bar_Vpeak;
        sigma_p_bar_VpeakM(i,j) = sigma_p_bar_Vpeak;
    end
end

%% Plots

%% Conductivity vs Work Function
for i = 1:length(Ncat_array)
    figure(100)
    semilogy(workfunction_LHS, sigma_n_bar_VpeakM(i, :))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
end  

for i = 1:length(Ncat_array)
    figure(101)
    semilogy(workfunction_LHS, sigma_p_bar_VpeakM(i, :))
    hold on
    xlabel('LHS workfunction [eV]')
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

workfunction_index = 9;
for i=1:241
    for j=1:length(Ncat_array)
    conductivity(j,i)=sigma_n_barM(j, workfunction_index,i);
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
for j = 1:length(workfunction_LHS)
    figure(201)
    semilogy(Vappt, squeeze(sigma_n_barM(3, j, :)))
    legstr_n2{j} = ['\Phi_l =', num2str(workfunction_LHS(j))];
    hold on
end

for j = 1:length(workfunction_LHS)
    figure(202)
    semilogy(Vappt, squeeze(sigma_p_barM(3, j, :)))
    legstr_p2{j} = ['\Phi_l =', num2str(workfunction_LHS(j))];
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

workfunction_index = 3;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    n_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3), 2);
    figure(203)
    semilogy(Vappt, n_int)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

for i = 1:length(Ncat_array)
    p_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3), 3);
    figure(204)
    semilogy(Vappt, p_int)
    legstr_p3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

figure(203)
xlabel('Voltage [V]')
ylabel('Electron density at the interface (cm-3)')
legend(legstr_n3)
hold off

figure(204)
xlabel('Voltage [V]')
ylabel('hole density interface (cm-3)')
legend(legstr_p3)
hold off

%% Plot electron and hole profiles at Vmax as a function of position

legstr_npx = {'', '', ''};
for i = 1:length(Ncat_array)
    dfplot.npx(sol_CV(i, workfunction_index), Vmax/k_scan);% Vmax/k_scan)
    legstr_npx{2*i-1 + 3} = ['n, Ncat =', num2str(Ncat_array(i))];
    legstr_npx{2*i + 3} = ['p, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_npx)
ylim([1e-1, 1e12])

%% Plot cation and anion densities at Vmax as a function of position




%dfplot.acx(sol_CV(3, 9), 3*(Vmax/k_scan)
%% Plot potential as a function position
workfunction_index = 9;
legstr_Vx = {'dielectric', 'interface', 'perovskite'};
for i = 1:length(Ncat_array)
    dfplot.Vx(sol_CV(i, workfunction_index), Vmax/k_scan);%Vmax/k_scan)
    legstr_Vx{i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_Vx)
%ylim([1e-1, 1e12])

%% Electron Modulability
idx = find(Vappt==0.4);
workfunction_index = 9;
legstr_n3 =[];
legstr_p3 =[];
 
for i = 1:length(Ncat_array)
    n_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3), 2);
    n_values_around_bp=n_int(idx-1:idx+1);
    Vappt_around_bp=Vappt(idx-1:idx+1);
    n_bp_modulability=gradient(n_values_around_bp,Vappt_around_bp);
    n_modulability_factor(:,i)= n_bp_modulability(2);
    n_store_log=log10(n_modulability_factor);
     figure(1111)
end
scatter(Ncat_array, n_store_log,'o', 'MarkerFaceColor', 'b');
set(gca,'xscale','log')

xlim([1e15 1e20])
ylim([4.6 8.5])
figure(1111)
xlabel('Cation concentration')
ylabel('Electron Modulability Factor (m_V_g)')
box on

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
%%

%% Plot ionic concentration at interface as function Vapp for different ion densities


%Plot similar for ions (instead of n_int put cat_int)
workfunction_index = 3;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    cat_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3)+1,4);
    logcat_int=log10(cat_int);
    figure(703)
    semilogy(Vappt,cat_int)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
figure(703)
xlabel('Voltage [V]')
ylabel('Cation Concentration (cm-3)')
legend(legstr_n3)
hold off
%% Modulability Ions

v_built_in=par.Phi_left-par.Phi_right;
workfunction_index = 3;
legstr_n3 =[];
legstr_p3 =[];
for i = 1:length(Ncat_array)
    cat_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3)+1,4);
    cat_int_log=log10(cat_int);
    figure(112)
    plot(Vappt, cat_int_log)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
figure(112)
xlabel('Voltage [V]')
ylabel('Cation Concentration (cm-3)')

hold off



%% Find how to get the contour done
x=workfunction_LHS;
y=Ncat_array;
z=sigma_n_bar_VpeakM;
z_log=log10(z);
surf(x,y,z_log);
%set(gca,'ZScale','log')
xlabel('Workfunction'), ylabel('Cation Concentration'), zlabel('Peak conductivity')
set(gca,'YScale','log')

%contour3();
%% Conductivity profiles
% So systematically you could look at the following.

% 1)Electrode workfunctions
% 2)Ion density 
% 3)1 ion, opposite charge (i.e. mobile anions)
% 4)2 ions
% 5)Different ion densities

% Ideally you would set some of these parameter explorations up as loops and extract peak conductivity then plot 
% on a contour plot with x = Ion density, y = Electrode workfunctions, z = peak conductivity for example.
