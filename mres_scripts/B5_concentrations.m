% Load device and edit sc_r value
par = pc('1_layer_MAPI_ITO_Ag.csv');

% CV parameters
V0 = 0; % starting voltage
Vmax = 1.2; % max voltage
Vmin = -1.2; % min voltage
cycles = 3; % number of CV cycles
points = (Vmax-Vmin)*200*cycles+1; % to get 100 points per volt

% equilibrate device
sol_eq = equilibrate(par);

% define parameter search arrays
sc_arr = logspace(-10,-6,5);    % surf rec vel
k_arr = logspace(-3,3,7);       % scan rate

% for ii = 1:length(sc_arr)
    par.sc_l = sc_arr(ii); % ITO electrode
    par.sc_r = sc_arr(ii); % Ag electrode
    
    par = refresh_device(par);

%       for jj = 1:length(k_arr)




title_txt = sprintf('%i cycles, sc_r = %g cms-1, sc_l = %g cms-1, rate = %g Vs-1'...
    ,[cycles, par.sc_r, par.sc_l, k_arr(jj)]);

sol_CV(ii) = doCV(


%       end
% end