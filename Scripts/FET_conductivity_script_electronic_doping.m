%% Line258
%try to see if yoou can plot for the time
%% Figure out the reason for the isobestic point
%Try for higher scan rate

%% Code pupose
% To compare the value of electron accumulation with ions and without ions and try to explain conductivity variance and some changes

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');

par = par_alox;     % Create temporary parameters object for overwriting parameters in loop

%% Initialise the parameter arrays
Ncat_array = logspace(16, 19, 4); %16, 19, 4
workfunction_LHS = -5.5:0.1:-4.2; %-5.5:0.1:-4.2;
%% Number of species
% par.N_ionic_species=2; %uncomment this to simulate with 2 ionic species
%% Set ionic mobility to be zero to simulate without affect of ions
% par.mu_c(:) = 0;
% par.mu_a(:) = 0;
%% For loop
for i = 1:length(Ncat_array)
    
    par.Ncat(:) = Ncat_array(i); %All layers get same ionic density. 
    par.Nani(:) = Ncat_array(i); %This is done for code purpose. 
   
    disp(['Cation density = ', num2str(Ncat_array(i)), ' cm^-3']);%num2str=Convert numbers to character representation
    for j = 1:length(workfunction_LHS) %loop to run for different electrode workfunction
        
        par.Phi_left = workfunction_LHS(j);
        disp(['LHS electrode workfunction = ', num2str(workfunction_LHS(j)), ' eV']);
        
        par = refresh_device(par);      % This line is required to rebuild various arrays used DF
        
        %% Find equilibrium
        soleq(i, j) = equilibrate(par);% Find and 
        dfplot.acx(soleq(i, j).ion) %plot the equilibirum solutions
        
        %% Current-voltage scan
        k_scan = 0.1;%Can check dependence of results on k_scan
        Vmax = 1.2;
        Vmin = -1.2;
%(from 0 to 1.2 to -1.2 to 0. Therefore (1.2x4/0.001)=(4800 scan points (checked in sol_CV(1, 1).t))
        tpoints=(20*(Vmax-Vmin)/k_scan)+1;
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        %tpoints is the No. of points in output time array
        sol_CV(i, j) = doCV(soleq(i, j).ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);%How is the number of time points determined?
       %Picked because it is a sensible divisor. Check for your voltage
       %range and then get it.
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

sigma_n_barM = zeros(length(Ncat_array), length(workfunction_LHS), length(sol_CV(1,1).t)); %Conductivity profiles for different Ncat 
sigma_p_barM = zeros(length(Ncat_array), length(workfunction_LHS), length(sol_CV(1,1).t)); %and Workfunctions for each different time steps
sigma_n_bar_VpeakM = zeros(length(Ncat_array), length(workfunction_LHS)); %Mean conducitvity across all the times
sigma_p_bar_VpeakM = zeros(length(Ncat_array), length(workfunction_LHS)); 

for i = 1:length(Ncat_array)
    for j = 1:length(workfunction_LHS)
       
        [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV(i,j));%call this function
        sigma_n_barM(i,j,:) = sigma_n_bar; %Put in calculated values for all times for different Ncat and wf
        sigma_p_barM(i,j,:) = sigma_p_bar; %A length(Ncat) x length(Wf) matrix for all time values will be created 
        %that has the mean conductivity
        sigma_n_bar_VpeakM(i,j) = sigma_n_bar_Vpeak; %This is calculated in sigma_ana
        sigma_p_bar_VpeakM(i,j) = sigma_p_bar_Vpeak;
        
    end
end

%% Plots
for i = 1:length(Ncat_array)
    figure(100)
    semilogy(workfunction_LHS, sigma_n_bar_VpeakM(i, :))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
    xlim([workfunction_LHS(1) workfunction_LHS(length(workfunction_LHS))])
    
end  

for i = 1:length(Ncat_array)
    figure(101)
    semilogy(workfunction_LHS, sigma_p_bar_VpeakM(i, :))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak hole conductivity [S cm-1]')
    legstr_p{i} = ['Ncat =', num2str(Ncat_array(i))];
    xlim([workfunction_LHS(1) workfunction_LHS(length(workfunction_LHS))])
    
end  
figure(100)
legend(legstr_n)
hold off
figure(101)
legend(legstr_p)
hold off


%% Plot average conductivity
for j = 1:length(workfunction_LHS)
    figure(201)
    semilogy(Vappt, squeeze(sigma_n_barM(3, j, :))) %remove singleton dimensions (reason?)
    legstr_n2{j} = ['\Phi_l =', num2str(workfunction_LHS(j))];
    xlim([Vmin Vmax]);
    hold on
end

for j = 1:length(workfunction_LHS)
    figure(202)
    semilogy(Vappt, squeeze(sigma_p_barM(3, j, :)))
    legstr_p2{j} = ['\Phi_l =', num2str(workfunction_LHS(j))];
    xlim([Vmin Vmax]);
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
workfunction_index =1;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    n_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3), 2); %n_int takes the electron density at the interface for all times, at space
    figure(203)
    semilogy(Vappt, n_int)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end

for i = 1:length(Ncat_array)
    p_int = sol_CV(i, workfunction_index).u(:, par.pcum0(3), 3);%par.pcum0(3) that corresponds to MAPI for variable number 2 that is the hole density. 

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
workfunction_index =1;
legstr_npx = {'', '', ''};
Vmax_t_index = find(sol_CV(1, 1).t == 80);
V_INDEX=Vmax/k_scan; %Better way to use the index?
for i = 1:length(Ncat_array)
    dfplot.npx(sol_CV(i, workfunction_index), Vmax/k_scan);% Vmax/k_scan) %this is the time index if I am correct
    legstr_npx{2*i-1 + 3} = ['n, Ncat =', num2str(Ncat_array(i))];
    legstr_npx{2*i + 3} = ['p, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_npx)
ylim([1e-1, 1e12])

%% Plot cation and anion densities at Vmax as a function of position
legstr_acx = {'', '', ''};
for i = 1:length(Ncat_array)
    dfplot.acx(sol_CV(i, workfunction_index), Vmax/k_scan)
    legstr_acx{2*i-1 + 3} = ['a, Ncat =', num2str(Ncat_array(i))];
    legstr_acx{2*i + 3} = ['c, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_acx)
%ylim([1e-1, 1e12])

%% Plot potential as a function position
workfunction_index =1;
legstr_Vx = {'dielectric', 'interface', 'perovskite'};
for i = 1:length(Ncat_array)
    dfplot.Vx(sol_CV(i, workfunction_index), 0);%Vmax/k_scan)
    legstr_Vx{i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_Vx)
%ylim([1e-1, 1e12])
%% Plot cation density as a function of voltage at the interface
cation_index=3;
legstr_n3 =[];
for i = 1:5
    cation_density_interface = sol_CV(cation_index, i).u(:, par.pcum0(3)+1,4); %Whats the value here for cation? Should par.pcum0(3) be different?
    figure(205)
    plot(Vappt, cation_density_interface)
    legstr_n3{i} = ['Phi =', num2str(workfunction_LHS(i))];
    hold on
end

figure(205)
xlabel('Voltage [V]')
ylabel('Cation density interface (cm-3)')
legend(legstr_n3)
hold off
%% Trying some dfplots
%   
for i = 4:-1:1
    dfplot.Vxacx(sol_CV(i,1), 0)
    subplot(2,1,1);
    hold on
    subplot(2,1,2)
    hold on
end
 subplot(2,1,1);
    hold off
    subplot(2,1,2)
    hold off
%Why does the graph disappear?

%% Potential due to ionic effect

   % dfplot.Vionxacx(sol_CV(1,7), 0)
   
   %% DFPLOT
  
   
dfplot.acx(sol_CV(3, 1),180);
%this is not time. check again
%% Plot individual values
%Ask how makemovie works
%  makemovie(sol_CV, @dfplot.npx, 0, [0, 1.5e18], 'npx', true, true);


%% Conductivity profiles
% So systematically you could look at the following.

% 1)Electrode workfunctions
% 2)Ion density 
% 3)1 ion, opposite charge (i.e. mobile anions)
% 4)2 ions
% 5)Different ion densities

% Ideally you would set some of these parameter explorations up as loops and extract peak conductivity then plot 
% on a contour plot with x = Ion density, y = Electrode workfunctions, z = peak conductivity for example.
%% Questions

%Why is dfana used for generating the voltage steps? Can it be done
%directly instead?
