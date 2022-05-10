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
Ncat_array = logspace(16, 19, 4);
workfunction_RHS = -4.5:0.1:-4.2;

%% while
for i = 1:length(Ncat_array)
    
    par.Ncat(:) = Ncat_array(i);
    par.Nani(:) = Ncat_array(i);
    
    disp(['Cation density = ', num2str(Ncat_array(i)), ' cm^-3']);
    for j = 1:length(workfunction_RHS) %loop to run for different electrode workfunction
        
        par.Phi_right = workfunction_RHS(j);
        disp(['RHS electrode workfunction = ', num2str(workfunction_RHS(j)), ' eV']);
        
        par.EF0(3) = workfunction_RHS(j);
        disp(['MAPI workfunction = ', num2str(workfunction_RHS(j)), ' eV']);
        
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
    end %loop runs till -5.5 to 4.9   
end

%% Analysis
Vappt = dfana.calcVapp(sol_CV(1,1));
% Preallocation
sigma_n_barM = zeros(length(Ncat_array), length(workfunction_RHS), length(sol_CV(1,1).t)); 
sigma_p_barM = zeros(length(Ncat_array), length(workfunction_RHS), length(sol_CV(1,1).t)); 
sigma_n_bar_VpeakM = zeros(length(Ncat_array), length(workfunction_RHS)); 
sigma_p_bar_VpeakM = zeros(length(Ncat_array), length(workfunction_RHS)); 

for i = 1:length(Ncat_array)
    for j = 1:length(workfunction_RHS)
        [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV(i,j));
        sigma_n_barM(i,j,:) = sigma_n_bar;
        sigma_p_barM(i,j,:) = sigma_p_bar;
        sigma_n_bar_VpeakM(i,j) = sigma_n_bar_Vpeak;
        sigma_p_bar_VpeakM(i,j) = sigma_p_bar_Vpeak;
    end
end

%% Plots
for i = 1:length(Ncat_array)
    figure(100)
    semilogy(workfunction_RHS, sigma_n_bar_VpeakM(i, :))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
end  

for i = 1:length(Ncat_array)
    figure(101)
    semilogy(workfunction_RHS, sigma_p_bar_VpeakM(i, :))
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


%% Plot average conductivity
for j = 1:length(workfunction_RHS)
    figure(201)
    semilogy(Vappt, squeeze(sigma_n_barM(3, j, :)))
    legstr_n2{j} = ['\Phi_l =', num2str(workfunction_RHS(j))];
    hold on
end

for j = 1:length(workfunction_RHS)
    figure(202)
    semilogy(Vappt, squeeze(sigma_p_barM(3, j, :)))
    legstr_p2{j} = ['\Phi_l =', num2str(workfunction_RHS(j))];
    hold on
end

figure(201)
xlabel('Voltage [V]')
ylabel('Average electron conductivity [Siemens]')
legend(legstr_n2)
hold off

figure(202)
xlabel('Voltage [V]')
ylabel('Average hole conductivity [Siemens]')
legend(legstr_p2)
hold off


%% Plot carrier concentration at interface as function Vapp for different ion densities
workfunction_index = 3;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    n_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3) + 1, 2);
    figure(203)
    semilogy(Vappt, n_int)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

for i = 1:length(Ncat_array)
    p_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3) + 1, 3);
    figure(204)
    semilogy(Vappt, p_int)
    legstr_p3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

figure(203)
xlabel('Voltage [V]')
ylabel('electron density interface (cm-3)')
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
    dfplot.npx(sol_CV(i, workfunction_index), round(3*Vmax/k_scan, 0));
    legstr_npx{2*i-1 + 3} = ['Ncat =', num2str(Ncat_array(i))];
    legstr_npx{2*i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_npx)
ylim([1e-1, 1e12])

%% Plot cation and anion densities at Vmax as a function of position
legstr_acx = {'', '', ''};
for i = 1:length(Ncat_array)
    dfplot.acx(sol_CV(i, workfunction_index), round(3*Vmax/k_scan, 0))
    legstr_acx{2*i-1 + 3} = ['Ncat =', num2str(Ncat_array(i))];
    legstr_acx{2*i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_acx)
%ylim([1e-1, 1e12])

%% Plot potential as a function position
legstr_Vx = {'dielectric', 'interface', 'perovskite'};
for i = 1:length(Ncat_array)
    dfplot.Vx(sol_CV(i, workfunction_index), 0);%Vmax/k_scan)
    legstr_Vx{i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_Vx)
%ylim([1e-1, 1e12])

%% Plot individual values

% makemovie(sol_CV, @dfplot.npx, 0, [0, 1.5e18], 'npx', true, true);

%% Conductivity profiles
% So systematically you could look at the following.

% 1)Electrode workfunctions
% 2)Ion density 
% 3)1 ion, opposite charge (i.e. mobile anions)
% 4)2 ions
% 5)Different ion densities

% Ideally you would set some of these parameter explorations up as loops and extract peak conductivity then plot 
% on a contour plot with x = Ion density, y = Electrode workfunctions, z = peak conductivity for example.
