% Load device and edit sc_r value
par = pc('1_layer_MAPI_ITO_Ag');
par.sc_r = 1e-8;
par.Phi_left = par.Phi_right;
par = refresh_device(par);

% initial solution
sol_eq = equilibrate(par);

% JV and function parameters
V0 = 0;
Vmax = 1.2;
Vmin = -1.2;
k = 1e-1;
cycles = 5;
points = 2401;
xpos = 0;

% run el only 
sol_el = doCV(sol_eq.el, 0, V0, Vmax, Vmin, k, 1, 241);

% run full CV
sol = doCV(sol_eq.ion, 0, V0, Vmax, Vmin, k, cycles, points);


% function to differentiate cycles in JV plot
% function JtotVappC(sol, xpos)
    % Same function as JtotVapp, but cycles are colour coded to differentiate
    % distinct cycles and forward/reverse direction

    % Obtain point position from x position, split into directions
    xmesh = sol.x;
    ppos = getpointpos(xpos, xmesh);

    % get JV arrays
    J = dfana.calcJ(sol);
    Vapp = dfana.calcVapp(sol);

    % get voltage gradient array
    gradV = diff(Vapp);
    gradV = [gradV(1), gradV];

    % find switch in direction in voltage
    dir_fwd = gradV > 0;
    dir_switch = diff(dir_fwd);

    % forward to reverse switch indices
    [~,f2r] = find(dir_switch == -1);
    % reverse to forward switch indices
    [~,r2f] = find(dir_switch == 1);

    % Complete potentially incomplete r2f array with last Vapp index
    % and define boolean to include last forward scan in plot
    if length(r2f) < length(f2r)
        r2f = [r2f,length(Vapp)];
        return_ticket = false;
    else
        return_ticket = true;
    end

    % initiate figure
    figure

    % set up legend entries
    leg_loop = ["Forward", "Reverse"];

    % choose colormap for plot lines
    c_map = jet(length(f2r));

    % initialize loop variables    
    leg_entries = [];
    Vstart = 1;
    
    % loop over number of cycles
    for i = 1:length(f2r)
        % iteratively split Vapp and Jtot into forward and reverse sections
        Vfwd = Vapp(Vstart:f2r(i));
        Vrev = Vapp(f2r(i):r2f(i));
        Jfwd = J.tot(Vstart:f2r(i),ppos);
        Jrev = J.tot(f2r(i):r2f(i),ppos);

        % set starting index for next loop
        Vstart = r2f(i);

        % create and append new legend entries
        leg_entries = [leg_entries,strcat(leg_loop, {' '}, num2str(i))];

        % plot forward and reverse scans for each cycle
        plot(Vfwd, Jfwd, '--', 'Color', c_map(i,:))
        hold on
        plot(Vrev, Jrev, ':', 'Color', c_map(i,:))

    end

    % include final forward scan to V0
    if return_ticket
        % get last forward JV values
        Vreturn = Vapp(r2f(end):end);
        Jreturn = J.tot(r2f(end):end,ppos);
        
        % Add plot to figure
        plot(Vreturn, Jreturn, 'k--')
        
        % append last legend entry to existing array
        leg_entries = [leg_entries, strcat("Forward", {' '}, ...
            num2str(length(f2r)+1))];
    end
    hold off

    % document graph
    xlabel('Applied Voltage, Vapp [V]');
    ylabel('Current Density, J [A cm^{-2}]');
    set(legend,'FontSize',14);
    set(legend,'EdgeColor',[1 1 1]);
    legend(leg_entries, 'location', 'northwest')

% end