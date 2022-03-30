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
%Ncat_array = logspace(16, 19, 4); %16, 19, 4
workfunction_LHS = -5.5:0.1:-4.2; %-5.5:0.1:-4.2;
%% Number of species
% par.N_ionic_species=2; %uncomment this to simulate with 2 ionic species
%% Set ionic mobility to be zero to simulate without affect of ions
% par.mu_c(:) = 0;
% par.mu_a(:) = 0;
%% For loop
    for j = 1:length(workfunction_LHS) %loop to run for different electrode workfunction
        
        par.Phi_left = workfunction_LHS(j);
        disp(['LHS electrode workfunction = ', num2str(workfunction_LHS(j)), ' eV']);
        
        par = refresh_device(par);      % This line is required to rebuild various arrays used DF
        
        %% Find equilibrium
        soleq(j) = equilibrate(par);% Find and 
        dfplot.acx(soleq(j).el) %plot the equilibirum solutions
        
        %% Current-voltage scan
        k_scan = 0.001;%Can check dependence of results on k_scan
        Vmax = 1.2;
        Vmin = -1.2;
%(from 0 to 1.2 to -1.2 to 0. Therefore (1.2x4/0.001)=(4800 scan points (checked in sol_CV(1, 1).t))
        
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        %tpoints is the No. of points in output time array
        sol_CV(j) = doCV(soleq(j).el, 0, 0, Vmax, Vmin, k_scan, 1, 241);%How is the number of time points determined?
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

%% Analysis
Vappt = dfana.calcVapp(sol_CV(1)); 

sigma_n_barM = zeros(length(workfunction_LHS), length(sol_CV(1).t)); %Conductivity profiles for different Ncat 
sigma_p_barM = zeros(length(workfunction_LHS), length(sol_CV(1).t)); %and Workfunctions for each different time steps
sigma_n_bar_VpeakM = zeros(length(workfunction_LHS)); %Mean conducitvity across all the times
sigma_p_bar_VpeakM = zeros(length(workfunction_LHS)); 


    for j = 1:length(workfunction_LHS)
       
        [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV(j));%call this function
        sigma_n_barM(j,:) = sigma_n_bar; %Put in calculated values for all times for different Ncat and wf
        sigma_p_barM(j,:) = sigma_p_bar; %A length(Ncat) x length(Wf) matrix for all time values will be created 
        %that has the mean conductivity
        sigma_n_bar_VpeakM(j) = sigma_n_bar_Vpeak; %This is calculated in sigma_ana
        sigma_p_bar_VpeakM(j) = sigma_p_bar_Vpeak;
        
    end


%% Plots
for j=1:length(workfunction_LHS)
    figure(500)
    semilogy(workfunction_LHS, sigma_n_bar_VpeakM(:,j))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak electron conductivity [S cm-1]')
     legend('Ncat =0');
    xlim([workfunction_LHS(1) workfunction_LHS(length(workfunction_LHS))])

end

for j=1:length(workfunction_LHS)
    figure(501)
    semilogy(workfunction_LHS, sigma_p_bar_VpeakM(:,j))
    hold on
    xlabel('LHS workfunction [eV]')
    ylabel('Peak hole conductivity [S cm-1]')
    legend('Ncat =0');
    xlim([workfunction_LHS(1) workfunction_LHS(length(workfunction_LHS))])
end


%% Plot average conductivity
for j = 1:length(workfunction_LHS)
    figure(502)
    semilogy(Vappt, sigma_n_barM(j,:)) %remove singleton dimensions (reason?)
  
    xlabel('Voltage [V]')
    ylabel('Average electron conductivity [Siemens]')
    legend('Ncat =0');
    xlim([Vmin Vmax]);
    hold on
end
hold off
for j = 1:length(workfunction_LHS)
    figure(503)
    semilogy(Vappt, sigma_p_barM(j,:))
    
     xlabel('Voltage [V]')
    ylabel('Average hole conductivity [Siemens]')
    legend('Ncat =0');
    xlim([Vmin Vmax]);
    xlim([Vmin Vmax]);
    hold on
end
hold off


%% Plot carrier concentration at interface as function Vapp for different ion densities
MAPI_index_and_rhs_wf_index=7;
cat_index=1;
    n_int = sol_CV(MAPI_index_and_rhs_wf_index).u(:, par.pcum0(3), 2);
    
    
    figure(504)
    semilogy(Vappt, n_int)
    legend('Ncat = 0')

    p_int = sol_CV(MAPI_index_and_rhs_wf_index).u(:, par.pcum0(3), 3);


    figure(505)
    semilogy(Vappt, p_int)
    legend('Ncat = 0')
   

figure(504)
xlabel('Voltage [V]')
ylabel('electron density interface (cm-3)')

hold off

figure(505)
xlabel('Voltage [V]')
ylabel('hole density interface (cm-3)')

hold off

%% Plot electron and hole profiles at Vmax as a function of position
workfunction_index =1;
legstr_npx = {'', '', ''};

V_INDEX=Vmax/k_scan; %Better way to use the index?

    dfplot.npx(sol_CV(workfunction_index), 0);% Vmax/k_scan) %this is the time index if I am correct
    
ylim([1e-1, 1e12])



%% Plot potential as a function position

workfunction_index =1;
legstr_Vx = {'dielectric', 'interface', 'perovskite'};

    dfplot.Vx(sol_CV( workfunction_index), 0);%Vmax/k_scan)
   
%ylim([1e-1, 1e12])

%% Trying some dfplots
 
% for i = 4:-1:1
%     dfplot.Vxacx(sol_CV(i,1), 0)
%     subplot(2,1,1);
%     hold on
%     subplot(2,1,2)
%     hold on
% end
%  subplot(2,1,1);
%     hold off
%     subplot(2,1,2)
%     hold off


%% Potential due to ionic effect

   % dfplot.Vionxacx(sol_CV(1,7), 0)
   

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
