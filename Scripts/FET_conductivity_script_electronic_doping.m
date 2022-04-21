%% Code pupose
% To check how doping affects things

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');

par = par_alox;     % Create temporary parameters object for overwriting parameters in loop

%% Initialise the parameter arrays
Ncat_array = logspace(16, 19, 4); %16, 19, 4
MAPI_Ef0 = -5.5:0.1:-4.2; %-5.5:0.1:-4.2;
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
    for j = 1:length(MAPI_Ef0) %loop to run for different electrode workfunction
        
        %par.Phi_right = MAPI_Ef0(j);%Change workfunction of right electrode with MAPI
        par.EF0(3)= MAPI_Ef0(j);%Change workfunction of MAPI
         %disp(['RHS electrode workfunction = ', num2str(MAPI_Ef0(j)), ' eV']);
         disp(['MAPI Ef0 = ', num2str(MAPI_Ef0(j)), ' eV']);
        
        par = refresh_device(par);      % This line is required to rebuild various arrays used DF
        
        %% Find equilibrium
        soleq(i, j) = equilibrate(par);% Find the equilibrium solutions 
        dfplot.acx(soleq(i, j).ion) %plot the equilibirum solutions
        
        %% Current-voltage scan
        k_scan = 0.001;
        Vmax = 1.2;
        Vmin = -1.2;
        tpoints=(2*(Vmax-Vmin)/(10*k_scan))+1;
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        %tpoints is the No. of points in output time array
        sol_CV(i, j) = doCV(soleq(i, j).ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
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
Vappt = dfana.calcVapp(sol_CV(1,1)); % Voltage applied on the device as a function of time
%If scan rate is 0.001 then why Vappt is not the same?
% Preallocation

sigma_n_barM = zeros(length(Ncat_array), length(MAPI_Ef0), length(sol_CV(1,1).t)); %Conductivity profiles for different Ncat 
sigma_p_barM = zeros(length(Ncat_array), length(MAPI_Ef0), length(sol_CV(1,1).t)); %and Workfunctions for each different time steps
sigma_n_bar_VpeakM = zeros(length(Ncat_array), length(MAPI_Ef0)); %Mean conducitvity acreoss all the times
sigma_p_bar_VpeakM = zeros(length(Ncat_array), length(MAPI_Ef0)); 

for i = 1:length(Ncat_array)
    for j = 1:length(MAPI_Ef0)
        [sigma_n_bar, sigma_p_bar, sigma_n_bar_Vpeak, sigma_p_bar_Vpeak] = sigma_ana(sol_CV(i,j));%call this function
        sigma_n_barM(i,j,:) = sigma_n_bar; %Put in calculated values for all times for different Ncat and wf
        sigma_p_barM(i,j,:) = sigma_p_bar; %A length(Ncat) x length(Wf) matrix for all time values will be created that has the mean conductivity
        sigma_n_bar_VpeakM(i,j) = sigma_n_bar_Vpeak; %This is calculated in sigma_ana
        sigma_p_bar_VpeakM(i,j) = sigma_p_bar_Vpeak;
    end
end
%What is the value in sigma_n_barM
%% Plots
for i = 1:length(Ncat_array)
    figure(100)
    semilogy(MAPI_Ef0, sigma_n_bar_VpeakM(i, :))
    hold on
    xlabel('MAPI workfunction [eV]')
    ylabel('Peak electron conductivity [S cm-1]')
    legstr_n{i} = ['Ncat =', num2str(Ncat_array(i))];
    xlim([MAPI_Ef0(1) MAPI_Ef0(length(MAPI_Ef0))])
end  

for i = 1:length(Ncat_array)
    figure(101)
    semilogy(MAPI_Ef0, sigma_p_bar_VpeakM(i, :))
    hold on
    xlabel('MAPI workfunction [eV]')
    ylabel('Peak hole conductivity [S cm-1]')
    legstr_p{i} = ['Ncat =', num2str(Ncat_array(i))];
    xlim([MAPI_Ef0(1) MAPI_Ef0(length(MAPI_Ef0))])
end  
figure(100)
legend(legstr_n)
hold off
figure(101)
legend(legstr_p)
hold off


%% Plot average conductivity
for j = 1:length(MAPI_Ef0)
    figure(201)
    semilogy(Vappt, squeeze(sigma_n_barM(3, j, :))) %remove singleton dimensions (reason?)
    legstr_n2{j} = ['\Phi_l =', num2str(MAPI_Ef0(j))];
    %xlim([Vmin Vmax]);
    hold on
end

for j = 1:length(MAPI_Ef0)
    figure(202)
    semilogy(Vappt, squeeze(sigma_p_barM(3, j, :)))
    legstr_p2{j} = ['\Phi_l =', num2str(MAPI_Ef0(j))];
    %xlim([Vmin Vmax]);
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
MAPI_index =7;
legstr_n3 =[];
legstr_p3 =[];

for i = 1:length(Ncat_array)
    n_int = sol_CV(i, MAPI_index).u(:, par.pcum0(3)+1, 2); %shouldn't this be 1? Confirm once
    figure(203)
    semilogy(Vappt, n_int)
    legstr_n3{i} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
%n_int takes the electron density at the interface for all times, at space
%par.pcum0(3) that corresponds to MAPI for variable number 2 that is the hole density. 
for i = 1:length(Ncat_array)
    p_int = sol_CV(i, MAPI_index).u(:, par.pcum0(3)+1, 3);%[time,space, variable]. 
%1. Electron density 2. Hole density 3. Cation density (where 1 or 2 mobile ionic
%carriers are stipulated)
%4. Anion density (where 2 mobile ionic carriers are stipulated)
%The spatial mesh x.
%The time mesh t.
%The parameters object par.


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
MAPI_index =7;
time_Vmin=round(3*Vmax/k_scan, 0);
time_Vmax=round(Vmax/k_scan, 0);
legstr_npx = {'', '', ''};
for i = 1:length(Ncat_array)
    dfplot.npx(sol_CV(i, MAPI_index),time_Vmax);% Vmax/k_scan corresponds to time.Vmax=1.2,k_scan=0.001
    %However for maximum negative voltage that should be at 3*Vmax/k_scan
% or at 3600 seconds
%We see there is a slight difference in the values.
    legstr_npx{2*i-1 + 3} = ['n, Ncat =', num2str(Ncat_array(i))];
    legstr_npx{2*i + 3} = ['p, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_npx)
ylim([1e-1, 1e12])

%% Plot cation and anion densities at Vmax as a function of position
legstr_acx = {'', '', ''};
for i = 1:length(Ncat_array)
    dfplot.acx(sol_CV(i, MAPI_index), Vmax/k_scan)
    legstr_acx{2*i-1 + 3} = ['a, Ncat =', num2str(Ncat_array(i))];
    legstr_acx{2*i + 3} = ['c, Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_acx)
%ylim([1e-1, 1e12])

%% Plot potential as a function position
MAPI_index =7;
legstr_Vx = {'dielectric', 'interface', 'perovskite'};
for i = 1:length(Ncat_array)
    dfplot.Vx(sol_CV(i, MAPI_index), time_Vmax);%Vmax/k_scan)
    legstr_Vx{i + 3} = ['Ncat =', num2str(Ncat_array(i))];
    hold on
end
legend(legstr_Vx)
%ylim([1e-1, 1e12])

%% Plot individual values

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
%How is tpoints determined?
