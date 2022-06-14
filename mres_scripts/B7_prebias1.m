% Freezing ions at different prebias voltages after slow CV

% load parameters
par = pc('1_layer_MAPI_ITO_Ag.csv');

% Set high sc
% sc_arr = [1e-8, 1e-5, 1];
sc_arr = 1;

% Set slow scan rate
% k_arr = [1e-3, 1e-10];
k_arr = 1e-10;

% CV parameters
V0 = 0;
Vmax = 5;   % High DeltaV
Vmin = -1;
cycles = 2;
points = (Vmax-Vmin)*200*cycles+1; % 100 points per Volt

% mobile anions on/off
for ai = 2
    par.N_ionic_species = ai;

    mob_anions = string(par.N_ionic_species == 2);

    for si = 1:length(sc_arr)
        % change sc
        par.sc_r = sc_arr; % right side
%         par.sc_l = sc_arr; % left side

        % par.mucat = par.mucat * 1e3;
        % par.muani = par.muani * 1e3;

        % par.mucat = 1;
        % par.muani = 1;

         par = refresh_device(par);

         % equilibrium solution
         sol_eq(ai) = equilibrate(par);

         % loop scan rates
         for ki = 1:length(k_arr)
             k = k_arr(ki);

            % Slow CV
            sol_CV_f(ai) = doCV(sol_eq(ai).ion, 0, V0, Vmax, Vmin, k, cycles, points);

            Vapp = dfana.calcVapp(sol_CV_f(ai));

            figure
            try
                dfplot.JtotVappC(sol_CV_f(ai),0)
            catch
                dfplot.JtotVapp(sol_CV_f(ai),0)
            end

            title(sprintf(['%i cycles, sc_r = %g cms-1, sc_l = %g cms-1,',...
                'rate = %g Vs-1, mobile anions: %s'],...
                cycles, sol_CV_f(ai).par.sc_r, sol_CV_f(ai).par.sc_l, k, mob_anions))
            plotbrowser

            % concentrations
            ACt(sol_CV_f(ai),cycles,0);
            title(sprintf(['%i cycles, sc_r = %g cms-1, sc_l = %g cms-1,',...
                'rate = %g Vs-1, mobile anions: %s'],...
                cycles, sol_CV_f(ai).par.sc_r, sol_CV_f(ai).par.sc_l, k, mob_anions))

            plotbrowser

            % acx
            [~,tarr] = dfana.splitsol(sol_CV_f(ai));
            tlast = tarr(end);
            figure
            dfplot.acx(sol_CV_f(ai), tlast);
            title(sprintf('ion distribution before failure at t = %1.5g s',tlast))

         end
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
