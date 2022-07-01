%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('./Input_files/alox.csv');
par_for_ions = par_alox;     % Create temporary parameters object for overwriting parameters in loop

%% Rough value of capacitance
A=1;
epsilon=par_for_ions.epp0*par_for_ions.epp(3)*par_for_ions.e;
d=par_for_ions.d(3);
Capacitance_rough=(A*epsilon)/d;
%% Set up parameters
 par_for_ions.Ncat(:)=1e16;

par_for_ions.Phi_left = -4.9;
par_for_ions.Phi_right = -4.9;
       

 %% Find equilibrium
soleq = equilibrate(par_for_ions);
       
 
        %% Current-voltage scan
k_scan= 0.01;
Vmax = 5;
Vmin = -5;
tpoints=400;
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_with_ions = doCV(soleq.ion, 0, 0, Vmax, Vmin, k_scan, 1, tpoints);
        
      
     


%% Capacitance_Manan Analysis
Vappt = dfana.calcVapp(sol_CV_with_ions); 
[Ctotal, Celectronic, Cionic] = capacitance_ana(sol_CV_with_ions,Vappt);%call this function
c12=Cionic;
 figure(7444)
 plot(Vappt,abs(Celectronic))
 legend('Capacitance')
xlabel('V applied')
ylabel('Capacitances (F/cm^2)')
 


%% Plots




%