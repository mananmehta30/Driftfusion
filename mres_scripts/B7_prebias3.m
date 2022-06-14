% Freezing ions at different prebias voltages after slow CV

% load parameters
par = pc('1_layer_MAPI_ITO_Ag.csv');
% par.Phi_right = par.Phi_left;
% par = refresh_device(par);

% Set high sc
sc = 1;

% Set slow scan rate
k = 1e-5;   % cations case
% k = 1e-8;   % anions case

% CV parameters
V0 = 0;
Vmax = 1;   % High DeltaV
Vmin = -1;   % cations case
% Vmin = -2;  % anions case
cycles = 1;
ppV = 100;  % 100 points per Volt
points = (Vmax-Vmin)*ppV*2*cycles+1; 

% Set no. of mobile ions
% 0 no mobile ions
% 1 cations only
% 2 cations and anions
ai = 1;

% loop mobile anions on/off
% for ai = 1:2
    par.N_ionic_species = ai;   % update parameters
    par.sc_r = sc;      % right side
%     par.sc_l = sc;    % left side
    par = refresh_device(par);

    % create true/false string for title etc.
    mob_anions = string(par.N_ionic_species == 2); 
    
    % equilibrium solution
    sol_eq = equilibrate(par);

%% Slow CV
sol_CV = doCV(sol_eq.ion, 0, V0, Vmax, Vmin, k, cycles, points);

% plot title text string
ti_txt = sprintf(['%i cycles, sc_r = %g cms-1, sc_l = %g cms-1,',...
'rate = %g Vs-1, mobile anions: %s'],...
cycles, sol_CV.par.sc_r, sol_CV.par.sc_l, k, mob_anions);

figure
try
dfplot.JtotVappC(sol_CV,0)
catch
dfplot.JtotVapp(sol_CV,0)
end

title(sprintf(['%i cycles, sc_r = %g cms-1, sc_l = %g cms-1,',...
'rate = %g Vs-1, mobile anions: %s'],...
cycles, sol_CV.par.sc_r, sol_CV.par.sc_l, k, mob_anions))
plotbrowser

% concentrations
ACt(sol_CV,cycles,0);
title(sprintf(['%i cycles, sc_r = %g cms-1, sc_l = %g cms-1,',...
'rate = %g Vs-1, mobile anions: %s'],...
cycles, sol_CV.par.sc_r, sol_CV.par.sc_l, k, mob_anions))

plotbrowser

%% Integration failure check

% Obtain final time/voltage points in solution
[~,tarr] = dfana.splitsol(sol_CV);
Vapp = dfana.calcVapp(sol_CV);

% Compare expected and output final time points
if length(tarr) < length(sol_CV.t)
    % Obtain failure time and voltage
    tfail = tarr(end);
%     tfail = tfail;
    Vfail = Vapp(end);
    
    fprintf('\nFailure at t = %g s, V = %g V \n', [tfail, Vfail])
%     % Plot ion densities over device at failure
%     figure
%     dfplot.acx(sol_CV, tfail)
% 
%     title(sprintf('Failure at t = %g s, V = %g V', [tfail, Vfail]))
%     
%     plotbrowser
end

%% makemovie
% makemovie(sol_CV(1), @dfplot.acx, 0, [0, 1.5e18], 'acx_equilibrium', 1, 0)

%% ELx and acx at certain voltages
V_arr = [-1, 0, 0.18 , 0.2, 0.22, .5, 1]; % no anions case
% V_arr = [-2, -1, -.5, 0, 0.18, 0.2, 0.22, 0.4, 0.6, 0.8, 2, 5]; % yes anions case
% V_arr = [-1, 0, 1];

for ii = 1:length(V_arr)
    if V_arr(ii) >= 0
        t_arr(ii) = V_arr(ii)/k;
    else
        t_arr(ii) = (abs(V_arr(ii))+2*Vmax)/k;
    end
end

for ti = 1:length(t_arr)
    V_string = sprintf(', Voltage = %.2g V',V_arr(ti));
    sti_txt = strcat(ti_txt, V_string);
    figure
    subplot(1,2,1)
    dfplot.acx(sol_CV, t_arr(ti))
%     set(gca, 'yscale', 'log')
    subplot(1,2,2)
    dfplot.ELx(sol_CV, t_arr(ti))
    sgtitle(sti_txt)
    
    plotbrowser
end

%% Freezing ions at set voltages
% V_fr_arr = [0.1,0.15,0.2,0.25,0.3,0.5,1,2,5];
% V_fr_arr = linspace(0.1,5,50);
% V_fr_arr = 0.1;
tic
% V_fr_arr = linspace(Vmin,Vmax,(Vmax-Vmin)*2+1);
V_fr_arr = [-2, -1, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4, 1, 2];

i_fail_pre = [];
i_fail_0V = [];
% loop over multiple freeze voltages
for Vi = 1:length(V_fr_arr)
    % Which voltage
    Vpre = V_fr_arr(Vi);

    % check if no prebias (0V)
    if V_fr_arr(Vi) == 0
        % if no prebias do not run doJV (doesn't work with 0V)
        % turn off ion mobility directly
        sol_eq.ion.par.mobseti = 0;
        fprintf('\nRunning CV with ions frozen at %.2g V \n', Vpre)
        sol_CV_fr(Vi) = doCV(sol_eq.el, 0, 0, 5, -5, 0.01, 1, 2001);
        % turn ion mobility back on
        sol_eq.ion.par.mobseti = 1;
        
    % all other voltages
    else    
        % jump to voltage where ions are to be frozen
        fprintf('\nApplying prebias of %.2g V \n', Vpre)
        sol_pre(Vi) = doJV(sol_eq.ion, 1e-10, 2*ppV*abs(Vpre)+1, 0, 1, 0, Vpre, 1);
        
        % Get applied voltage array of JV forward scan
        Vapp = dfana.calcVapp(sol_pre(Vi).dk.f);
        Vfinal1 = Vapp(end);
    
        
        if Vfinal1 == Vpre
            % return to 0V with ion mobility turned off
            fprintf('\nReturnin to 0V from %.2g V \n', V_fr_arr(Vi))
            sol_0V(Vi) = doJV(sol_pre(Vi).dk.f, 0.1, 101, 0, 0, Vpre, 0, 1);
            
            Vapp = dfana.calcVapp(sol_0V(Vi).dk.f);
            Vfinal2 = Vapp(end);

            if Vfinal2 == 0
                %  run a CV with frozen ions
                fprintf('\nRunning CV with ions frozen at %.2g V \n', V_fr_arr(Vi))
                sol_CV_fr(Vi) = doCV(sol_0V(Vi).dk.f, 0, 0, 5, -5, 0.01, 1, 2001);
            else 
                fprintf('\nFailed to return to 0 V from prebias of %g V\n',Vpre)
                i_fail_0V = [i_fail_0V,Vi];
            end
            
        else % if forward JV didn't reach desired prebias
            fprintf('\nFailed to reach prebias of %g V\n',Vpre)
            i_fail_pre = [i_fail_pre,Vi];
            
        end
    end
end

%     % check if ions are frozen with 
%     t_arr = [0, 5, 10, 15, 20].*(1/0.01);
%     figure(100+ai)
%     subplot(3,3,Vi)
%     try
%         dfplot.acx(sol_CV_fr(Vi), t_arr)
%         title(sprintf('frozen at %.2g V',Vp))
%     end
toc
%%
for ei = 1:length(V_fr_arr)
    succ(ei) = any(find(~structfun(@isempty,sol_CV_fr(ei))));
end

V_succ = V_fr_arr(succ ~= 0);

%% plotting all freeze voltages
figure()
succ= [];
for Vi = 1:length(V_fr_arr)
    % plot CV
    hold on
    try
        dfplot.JtotVapp(sol_CV_fr(Vi), 0)
        succ = [succ,Vi];
    catch
        warning('\nNo plot output at %g V \n', V_fr_arr(Vi));
    end
    
end

% Record all successful freeze voltages
% V_succ = V_fr_arr(succ);

hold off

leg_txt = string(V_fr_arr(succ)) + " V";
legend(leg_txt)
title(sprintf('Frozen ions at voltages, mobile anions %s',mob_anions));
plotbrowser

%%
for si = 1:length(succ)
    % Only loop over successful attempts
    fi = succ(si);
    % Calculate resistance around 0V
    % Get Vapp and Jtot
    Vapp = dfana.calcVapp(sol_CV_fr(fi));
    J_arr = dfana.calcJ(sol_CV_fr(fi));
    Jtot = J_arr.tot(:,1);
    % Get index of zero crossing
    indices = find(~any(Vapp,1));
    ind(si) = indices(1);
    % Get JV gradient
    JV_grad = gradient(Jtot,Vapp);
    JV_grad0(si) = JV_grad(ind(si));


end

%%
figure
plot(V_succ, abs(JV_grad0), 'x')
set(gca,'yscale', 'log')
title('Conductance at 0 V after applied pre-bias')
xlabel('Pre-bias [V]')
ylabel('Conductance [\Omega^{-1}]')




% plotbrowser

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

% plot resistance and ion distribution as function of prebias voltage

% change workfunctions