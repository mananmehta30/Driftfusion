% Load device and edit sc_r value
par = pc('1_layer_MAPI_ITO_Ag');
par.sc_r = 1e-8;
par = refresh_device(par);
% initial solution
sol_eq = equilibrate(par);

% CV and function parameters
V0 = 0;
Vmax = 1.2;
Vmin = -1.2;
k = 1e-3;
cycles = 5;
points = 2400;
xpos = 0;

% run el only 
sol_el = doCV(sol_eq.el, 0, V0, Vmax, Vmin, k, 1, 241);

% run full CV
sol = doCV(sol_eq.ion, 0, V0, Vmax, Vmin, k, cycles, points);


%%
% function JtotVappC(sol, xpos, Vmax, Vmin, cycles, points)
            % Obtain point position from x position, split into directions
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);
            
            % points per cycle
            ppc = round((points - 1)/cycles);
            
            % dV
            dV = Vmax - Vmin;
            
            % steps in each part of each cycle
            sharepos = abs(Vmax/(Vmax-Vmin)/2);
            shareneg = abs(Vmin/(Vmax-Vmin)/2);
            
            cyclepos = round(sharepos * ppc);
            cycleneg = round(shareneg * ppc);
            
            % Set up for loop
            c_map = parula(cycles);
            
            % initialize legend entries
            leg_entries = [];
            loop_entries = ["fwd pos","rev pos","rev neg","fwd pos"];
            
            counter = 0;
            for i = 1:cycles
                % Define forward/reverse ranges
                fpos = (counter+1:counter+1+cyclepos);
                rpos = (counter+1+cyclepos:counter+2*cyclepos);
                rneg = (counter+2*cyclepos: counter+2*cyclepos+cycleneg);
                fneg = (counter+2*cyclepos+cycleneg:ppc*i);
                
                % Set counter for next loop
                counter = ppc*i;
                
                % plot all 4 directions for cycle i
                plot(Vapp(fpos), J.tot(fpos,ppos),'-', 'Color', c_map(i,:))
                hold on
                plot(Vapp(rpos), J.tot(rpos,ppos), '--', 'Color', c_map(i,:))
                plot(Vapp(rneg), J.tot(rneg,ppos), '-.', 'Color', c_map(i,:))
                plot(Vapp(fneg), J.tot(fneg,ppos), ':', 'Color', c_map(i,:))
                
                % create new legend entries
                cycle_no = sprintf('Cycle %i,',i);
                cycle_leg = strcat(cycle_no,{' '}, loop_entries);
                leg_entries = [leg_entries, cycle_leg];
                
                
            end
            
            % document graph
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend,'FontSize',14);
            set(legend,'EdgeColor',[1 1 1]);
            legend(leg_entries)
            
            hold off
% end
