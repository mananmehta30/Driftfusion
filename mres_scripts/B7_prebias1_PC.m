% Freezing ions at different prebias voltages after slow CV

% load parameters
par = pc('1_layer_MAPI_ITO_Ag.csv');

% Set high sc
sc_arr = [1e-10,1];

% Set slow scan rate
k_arr = [1e-3,1e-5,1e-10];

% CV parameters
V0 = 0;
Vmax = 5;   % High DeltaV
Vmin = -5;
cycles = 1;
ppV = 20;
points = (Vmax-Vmin)*ppV*cycles+1; % 100 points per Volt

for si = 1:length(sc_arr)
    % change sc    
    par.sc_r = sc_arr(si); % right side
    % par.sc_l = sc; % left side
    
    par.mucat = par.mue/1e3;
    % par.muani = par.muani * 1e3;

    % par.mucat = 1;
    % par.muani = 1;

     par = refresh_device(par);
     
     % equilibrium solution
     sol_eq(si) = equilibrate(par);
     
     % Reduce solver tolerance
     sol_eq(si).ion.par.AbsTol = 1e-4;
     sol_eq(si).ion.par.RelTol = 1e-10;
     
     % loop scan rates
     for ki = 1:length(k_arr)
         k = k_arr(ki);
         
        % Slow CV
        % sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
        sol_CV_f(si,ki) = doCV(sol_eq(si).ion, 0, V0, Vmax, Vmin, k, cycles, points);
        
        Vapp = dfana.calcVapp(sol_CV_f(si,ki));

        figure
%         try
%             dfplot.JtotVappC(sol_CV_f(si,ki),0)
%         catch
            dfplot.JtotVapp(sol_CV_f(si,ki),0)
%         end

        title(sprintf(['%i cycles, sc_r = %.0e cms-1, sc_l = %.0e cms-1, ',...
            'rate = %.0e Vs-1'], [cycles, par.sc_r, par.sc_l, k]))
        plotbrowser
        
        % concentrations
%         ACt(sol_CV_f(si,ki),cycles,0);
%         title(sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1, V_{max} = %.2f V, V_{min} = %.2f V',...
%         [cycles, sol_CV_f(si,ki).par.sc_r, sol_CV_f(si,ki).par.sc_l, k, Vmax, Vmin]))

     end
end




% CV forward & backward
% compare resulting curves
% should match to indicate equilibrium of ions

% Select voltage points

% pick out solution at voltages

% set ion mobility to 0

% return to 0 V (Will)

% Run CV with frozen ions

% measure IV slope around 0 V

% classify state (resistive/rectifying)

% plot resistanc and ion distribution as function of prebias voltage

% change workfunctions