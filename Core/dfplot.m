classdef dfplot
    % DRIFTFUSION Plotting class - contains methods for plotting
    %
    % List of available plots:
    % DFPLOT.JT = Currents as a function of time
    % DFPLOT.JX = Total currents as a function of position
    % DFPLOT.jx = Carrier fluxes as a function of position
    % DFPLOT.JV = Current-voltage curve using a solution from DOJV-
    % now outdated - preferrable to use DOCV and DFPLOT.JTOTVAPP
    % DFPLOT.JDDX = Drift and diffusion currents as a function of position
    % DFPLOT.VOCT = Open circuit voltage as a function of time
    % DFPLOT.PLT = Integrated radiative recombination rate as a function of
    % time
    % DFPLOT.VAPPT = Applied voltage as a function of time
    % DFPLOT.JVAPP = Current components as a function of the applied voltage at
    % position defined by XPOS
    % DFPLOT.JTOTVAPP = Total current as a function of the applied voltage at
    % position defined by XPOS
    % DFPLOT.LOGJVAPP = Current components as a function of the applied voltage
    % using log y-axis at position defined by XPOS
    % DFPLOT.XMESH = Plots the xmesh position vs point number
    % DFPLOT.VX = Electrostatic potential as a function of position
    % DFPLOT.NPX = Electron and hole densities as a function of position
    % DFPLOT.ACX = Anion and cation densities as a function of position
    % DFPLOT.GX = Generation rate as a function of position
    % DFPLOT.GXT = Generation rate as a function of position and time
    % DFPLOT.RX = Recombination rate components as a function of position
    % DFPLOT.JRECVAPP = Recomonbination components as a function of the applied
    % voltage
    % DFPLOT.FT = Electric field as a function of time
    % DFPLOT.SIGMAT = Integrated charge density in cm-2 as a function of time
    % DFPLOT.QT = Charge density in Coulombs cm-2 integrated between within the
    % range [X1, X2] as a function of time
    % DFPLOT.QVAPP = Charge density in Coulombs cm-2 integrated between within the
    % range [X1, X2] as a function of applied voltage
    % DFPLOT.RHOX = Volumetric charge density as a function of position
    % DFPLOT.DELTARHOX = Change in volumetric charge density as a function of position
    % DFPLOT.RHOXFXVX = Volumetric charge density, Electric field and Electrostatic potential
    % as a function of position- stacked plot
    % DFPLOT.RHOXVX = Volumetric charge density and Electrostatic potential as a function
    % of position- stacked plot
    % DFPLOT.ELX = Energy level diagram as a function of position
    % DFPLOT.ELXNPX = Energy level diagram, electron and hole densities
    % DFPLOT.ELXNPXACX = Energy level diagram, electron and hole densities,
    % and anion and cation densities, 3 panel, stacked
    % DFPLOT.VXACX = Electrostatic potential and anion and cation densities, 2
    % panel, stacked
    % DFPLOT.VIONXACX = Electrostatic potential due to ionic charge and anion and cation densities, 2
    % panel, stacked
    % DFPLOT.FIONT = Electric field due to the ionic charge as a function of
    % time

    % Plotting functions that are a function of position can accept a time
    % array as the second argument- the procedure will loop and plot the
    % solution at multiple times.
    % The third optional argument defines the x-range.
    % For plotting functions that are a function of time, the second argument
    % is generally the position at which the value is taken- see the comments
    % of individual methods below for further details

    %% LICENSE
    % Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
    % Imperial College London
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published
    % by the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    %% Start code
    methods (Static)

        function Jt(sol, xpos)
            % Currents as a function of time
            % SOL = solution structure
            % XPOS = the readout position
            [~,t,~,~,~,~,~,~,~,~] = dfana.splitsol(sol);

            [J, j, xmesh] = dfana.calcJ(sol);
            ppos = getpointpos(xpos, xmesh);

            
            plot(t, J.n(:, ppos),t, J.p(:, ppos),t, J.a(:, ppos),t, J.c(:, ppos), t, J.disp(:,ppos), t, J.tot(:, ppos));
            legend('Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtotal')
            xlabel('time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end

        function Jx(varargin)
            % Plots the current components
            % VARARGIN = [SOL, TARR, XRANGE]
            % SOL = Solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [XMIN, XMAX]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [J, j, x] = dfana.calcJ(sol);

            
            dfplot.x2d(sol, x, {J.n, J.p, J.a, J.c, J.disp, J.tot},...
                {'Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtot'}, {'-','-','-','-','-','-'},...
                'Current density [Acm-2]', tarr, xrange, 0, 0);
        end

        function jx(varargin)
            % Plots the carrier fluxes
            % VARARGIN = [SOL, TARR, XRANGE]
            % SOL = Solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [XMIN, XMAX]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [J, j, x] = dfana.calcJ(sol);

            dfplot.x2d(sol, par.x_ihalf, {j.n, j.p, j.a, j.c, j.disp},{'jn', 'jp', 'ja', 'jc', 'jdisp'},...
                {'-','-','-','-','-'}, 'Flux [cm-2 s-1]', tarr, xrange, 0, 0);
        end       
     
        function JV(JV, option)
            % JV - a solution from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs

            if option == 1 || option == 3
                J.dk.f = dfana.calcJ(JV.dk.f);
                Vapp.dk.f = dfana.calcVapp(JV.dk.f);
                J.dk.r = dfana.calcJ(JV.dk.r);
                Vapp.dk.r = dfana.calcVapp(JV.dk.r);

                
                plot(Vapp.dk.f, J.dk.f.tot(:,end), '--', Vapp.dk.r, J.dk.r.tot(:,end));
                hold on
            end

            if option == 2 || option == 3

                J.ill.f = dfana.calcJ(JV.ill.f);
                Vapp.ill.f = dfana.calcVapp(JV.ill.f);
                J.ill.r = dfana.calcJ(JV.ill.r);
                Vapp.ill.r = dfana.calcVapp(JV.ill.r);

                
                plot(Vapp.ill.f, J.ill.f.tot(:,end),'--')%, 'Color', [0, 0.4470, 0.7410]);
                hold on
                plot(Vapp.ill.r, J.ill.r.tot(:,end));%,'Color', [0, 0.4470, 0.7410]);
            end

            
            %ylim([-30e-3, 10e-3]);
            xlabel('Applied voltage [V]')
            ylabel('Current density [Acm-2]');
            hold off
        end
        
        function Jddx(varargin)
            % 
            % drift and diffusion currents as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [~, Jdd, x] = dfana.Jddxt(sol);
            
            
            dfplot.x2d(sol, x, {Jdd.ndiff, Jdd.ndrift, Jdd.pdiff, Jdd.pdrift,...
                Jdd.adiff, Jdd.adrift, Jdd.cdiff, Jdd.cdrift},...
                {'Jn,diff', 'Jn,drift', 'Jp,diff', 'Jp,drift', 'Ja,diff', 'Ja,drift', 'Jc,diff', 'Jc,drift'},...
                {'.','-','.','-','.','-','.','-'},'Current density [Acm-2]', tarr, xrange, 0, 0);
        end
        
        function Voct(sol)
            [~,t,~,~,~,~,~,~,~,~] = dfana.splitsol(sol);
            Voc = dfana.calcVQFL(sol);
            
            plot(t, Voc)
            xlabel('Time [s]')
            ylabel('Voc [V]')
        end
        
        function PLt(sol)
            [~,t,~,~,~,~,~,~,~,~] = dfana.splitsol(sol);
            PL = dfana.PLt(sol);
            
            plot(t, PL)
            xlabel('Time [s]')
            ylabel('PL [cm-2s-1]')
        end
        
        function Vappt(sol)
            [~,t,~,~,~,~,~,~,~,~] = dfana.splitsol(sol);
            % Difference in potential between the left and right boundary
            Vapp = dfana.calcVapp(sol);
            
            
            plot(t, Vapp);
            xlabel('Time [s]')
            ylabel('Vapp [V]')
        end
        
        function JVapp(sol, xpos)
            % Obtain point position from x position
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);

            
            plot(Vapp, J.n(:, ppos),Vapp, J.p(:, ppos), Vapp, J.c(:, ppos), Vapp, J.a(:, ppos),...
                Vapp, J.disp(:,ppos), '--', Vapp, J.tot(:, ppos), '--');
            legend('Jn', 'Jp', 'Jc', 'Ja', 'Jdisp', 'Jtotal')
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end

        function JtotVapp(sol, xpos)
            % Obtain point position from x position
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);

            
            semilogy(Vapp, abs(J.tot(:, ppos)));
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end

        function JtotVappC(sol, xpos)
            % Same function as JtotVapp, but cycles are colour coded to 
            % differentiate distinct cycles and forward/reverse direction

            % Obtain point position from x position
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

        end
%         Old cycle labelling JV plot function
%         function JtotVappC(sol, xpos, Vmax, Vmin, cycles, points)
%             % Obtain point position from x position, split into directions
%             % DOESNT WORK WELL WITH OTHER PLOTS IN SAME FIGURE
%             xmesh = sol.x;
%             ppos = getpointpos(xpos, xmesh);
% 
%             J = dfana.calcJ(sol);
%             Vapp = dfana.calcVapp(sol);
%             
%             % points per cycle
%             ppc = round((points - 1)/cycles);
%             
%             % dV
%             dV = Vmax - Vmin;
%             
%             % steps in each part of each cycle
%             sharepos = abs(Vmax/(Vmax-Vmin)/2);
%             shareneg = abs(Vmin/(Vmax-Vmin)/2);
%             
%             cyclepos = round(sharepos * ppc);
%             cycleneg = round(shareneg * ppc);
%             
%             % Set up for loop
%             c_map = parula(cycles);
%             
%             % initialize legend entries
%             leg_entries = [];
%             loop_entries = ["fwd pos","rev pos","rev neg","fwd pos"];
%             
%             counter = 0;
%             
%             % loop each cycle of the JV
%             for i = 1:cycles
%                 % Define forward/reverse ranges
%                 fpos = (counter+1:counter+1+cyclepos);
%                 rpos = (counter+1+cyclepos:counter+2*cyclepos);
%                 rneg = (counter+2*cyclepos: counter+2*cyclepos+cycleneg);
%                 fneg = (counter+2*cyclepos+cycleneg:ppc*i);
%                 
%                 % Set counter for next loop
%                 counter = ppc*i;
%                 
%                 % plot all 4 directions for cycle i
%                 plot(Vapp(fpos), J.tot(fpos,ppos),'-', 'Color', c_map(i,:))
%                 hold on
%                 plot(Vapp(rpos), J.tot(rpos,ppos), '--', 'Color', c_map(i,:))
%                 plot(Vapp(rneg), J.tot(rneg,ppos), '-.', 'Color', c_map(i,:))
%                 plot(Vapp(fneg), J.tot(fneg,ppos), ':', 'Color', c_map(i,:))
%                 
%                 % create new legend entries
%                 cycle_no = sprintf('Cycle %i,',i);
%                 cycle_leg = strcat(cycle_no,{' '}, loop_entries);
%                 leg_entries = [leg_entries, cycle_leg];
%                 
%                 
%             end
%             
%             % document graph
%             xlabel('Applied Voltage, Vapp [V]');
%             ylabel('Current Density, J [A cm^{-2}]');
%             set(legend,'FontSize',14);
%             set(legend,'EdgeColor',[1 1 1]);
%             legend(leg_entries,'location', 'eastoutside')
%             
%             hold off
%         end

        function logJVapp(sol, xpos)
            % plot the log of the mod J

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);

            
            semilogy(Vapp, abs(J.tot(:,ppos)), Vapp, abs(J.n(:,ppos)), Vapp, abs(J.p(:,ppos)), Vapp, abs(J.a(:,ppos)),Vapp, abs(J.c(:,ppos)), Vapp, abs(J.disp(:,ppos)));
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');
            legend('Jtot', 'Jn', 'Jp', 'Ja', 'Jc', 'Jdisp')
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end

        function logJVapp3D(sol, xpos, ylogon)

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            t = sol.t;
            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol)';
            Jtot=J.tot(:, ppos);

            
            surface('XData', [Vapp Vapp],             ... % N.B.  XYZC Data must have at least 2 cols
                'YData', [abs(Jtot) abs(Jtot)],             ...
                'ZData', [t' t'], ...
                'CData', [t' t'],             ...
                'FaceColor', 'none',        ...
                'EdgeColor', 'interp',      ...
                'Marker', 'none','LineWidth',1);
            s1 = gca;
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');

            if ylogon
                set(s1,'YScale','log');
            else
                set(s1,'YScale','linear');
            end
            hold off
        end

        function xmesh(sol)
            plot(sol.x)
            ylabel('Point')
            xlabel('Position [cm]')
        end

        function Vx(varargin)
            % Electrostatic potential as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            
            dfplot.x2d(sol, x, {V},{'V'},{'-'},'Electrostatic potential [V]', tarr, xrange, 0, 0);
        end

        function Fx(varargin)
            % Electrostatic potential as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            F = dfana.calcF(sol);
            
            
            dfplot.x2d(sol, x, {F},{'F'},{'-'},'Electric field [Vcm-1]', tarr, xrange, 0, 0);
        end
        
        function npx(varargin)
            % Carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-','-'},...
                'Carrier density [cm-3]', tarr, xrange, 0, 1)
        end

        function acx(varargin)
            % Ionic carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            dfplot.x2d(sol, x, {a,c},{'anion','cation'}, {'-','-'},...
                'Ionic carrier density [cm-3]', tarr, xrange, 0, 0);
        end

        function gx(varargin)
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [g1, g2, g] = dfana.calcg(sol);

            
            dfplot.x2d(sol, par.x_ihalf, {g1, g2, g}, {'g1', 'g2', 'g total'},...
                {'-','-','-'}, 'Generation rate [cm-3s-1]', tarr, xrange, 0, 0);
        end

        function gxt(sol)
            % Carrier densities as a function of position
            par = sol.par;
            [~,t,~,~,~,~,~,~,~,~] = dfana.splitsol(sol);
            [~, ~, g] = dfana.calcg(sol);
            xnm = par.x_ihalf*1e7;

            
            surf(xnm, t, g)
            xlabel('Position [cm]')
            ylabel('Time [s]')
            zlabel('Generation rate [cm^{-3}s^{-1}]')
        end

        function rx(varargin)
            % Recombination rates as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            r = dfana.calcr(sol);

            
            dfplot.x2d(sol, x, {r.btb, r.srh, r.tot},{'rbtb', 'rsrh', 'rtot'},...
                {'-','-','-'}, 'Recombination rate [cm-3s-1]', tarr, xrange, 0, 0);
        end

        function rsrhx(varargin)
            % Recombination rates as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            r = dfana.calcr(sol);

            
            dfplot.x2d(sol, x, {r.srh},{''},...
                {'-'}, 'SRH recombination rate [cm-3s-1]', tarr, xrange, 0, 1);
        end
        
        function JrecVapp(JV, option)
            % Plots recombination currents for JV

            % JV - a solution from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs

            if option == 1 || option == 3
                J.dk.f = dfana.calcJ(JV.dk.f);
                Vapp.dk.f = dfana.calcVapp(JV.dk.f);
                J.dk.r = dfana.calcJ(JV.dk.r);
                Vapp.dk.r = dfana.calcVapp(JV.dk.r);

                
                plot(Vapp.dk.f, J.dk.f.tot(:,end), '--', Vapp.dk.r, J.dk.r.tot(:,end));
                hold on
            end

            if option == 2 || option == 3
                solf = JV.ill.f;
                solr = JV.ill.r;
                par = solf.par;
                pcum0 = par.pcum0;

                J.ill.f = dfana.calcJ(JV.ill.f);
                Vapp.ill.f = dfana.calcVapp(JV.ill.f);
                J.ill.r = dfana.calcJ(JV.ill.r);
                Vapp.ill.r = dfana.calcVapp(JV.ill.r);

                r_f = dfana.calcr(JV.ill.f);
                Jrec_btb_f = JV.ill.f.par.e*trapz(JV.ill.f.x, r_f.btb, 2);
                Jrec_srhint_f = JV.ill.f.par.e*trapz(JV.ill.f.x(pcum0(2)+1:pcum0(3)), r_f.srh(:,pcum0(2)+1:pcum0(3)), 2)...
                    +JV.ill.f.par.e*trapz(JV.ill.f.x(pcum0(4)+1:pcum0(5)), r_f.srh(:,pcum0(4)+1:pcum0(5)), 2);
                Jrec_srhbulk_f = JV.ill.f.par.e*trapz(JV.ill.f.x(pcum0(3)+1:pcum0(4)), r_f.srh(:,pcum0(3)+1:pcum0(4)), 2);
                Jrec_tot_f = JV.ill.f.par.e*trapz(JV.ill.f.x, r_f.tot, 2);

                r_rev = dfana.calcr(JV.ill.r);
                Jrec_btb_r = JV.ill.f.par.e*trapz(JV.ill.r.x, r_rev.btb, 2);
                Jrec_srhint_r = JV.ill.r.par.e*trapz(JV.ill.r.x(pcum0(2)+1:pcum0(3)), r_rev.srh(:,pcum0(2)+1:pcum0(3)), 2)...
                    +JV.ill.r.par.e*trapz(JV.ill.r.x(pcum0(4)+1:pcum0(5)), r_rev.srh(:,pcum0(4)+1:pcum0(5)), 2);
                Jrec_srhbulk_r = JV.ill.r.par.e*trapz(JV.ill.r.x(pcum0(3)+1:pcum0(4)), r_rev.srh(:,pcum0(3)+1:pcum0(4)), 2);
                Jrec_tot_r = JV.ill.r.par.e*trapz(JV.ill.r.x, r_rev.tot, 2);

                cc=lines(4);

                
                plot(Vapp.ill.f, J.ill.f.tot(:,end),'--', 'Color', cc(1,:));
                hold on
                plot(Vapp.ill.r, J.ill.r.tot(:,end), 'Color', cc(1,:));
                % Recombination currents
                plot(Vapp.ill.f, Jrec_btb_f,'--','Color', cc(2,:));
                plot(Vapp.ill.r, Jrec_btb_r,'Color', cc(2,:));
                plot(Vapp.ill.f, Jrec_srhint_f,'--','Color', cc(3,:));
                plot(Vapp.ill.r, Jrec_srhint_r, 'Color', cc(3,:));
                plot(Vapp.ill.f, Jrec_srhbulk_f,'--','Color', cc(4,:));
                plot(Vapp.ill.r, Jrec_srhbulk_r, 'Color', cc(4,:));
            end

            
            ylim([-30e-3, 10e-3]);
            xlabel('Applied voltage [V]')
            ylabel('Current density [Acm-2]');
            hold off
            legend('Illumated for', 'Illumated rev','Jrec,btb for','Jrec,btb rev'...
                ,'Jrec,srh-int for','Jrec,srh-int rev','Jrec,srh-bulk for','Jrec,srh-bulk rev')
        end

        function Ft(sol, xpos)
            % Absolute field strength F as a function of time at point
            % position XPOS
            [~,t,xmesh,~,~,~,~,~,~,~] = dfana.splitsol(sol);
            ppos = getpointpos(xpos, xmesh);

            F = dfana.calcF(sol);

            
            plot(t, F(:,ppos))
            xlabel('Time [s]')
            ylabel(['Electric Field at pos x = ', num2str(round(xpos*1e7)), 'nm [Vcm-1]'])
        end

        function sigmat(sol)
            % Plot the integrated space charge density [cm-2] as a function of time
            sigma = dfana.calcsigma(sol);
            [~,t,~,~,~,~,~,~,~,~] = dfana.splitsol(sol);  
            
            plot(t, sigma)
            xlabel('Time [s]')
            ylabel('sigma [C cm-2]')
        end

        function Qt(sol, x1, x2)
            % Plot the integrated space charge density in Coulombs [Ccm-2] as a function of time
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            p1 = find(x<=x1);
            p1 = p1(end);
            p2 = find(x<=x2);
            p2 = p2(end);

            rho = dfana.calcrho(sol);
            Q = par.e*trapz(x(p1:p2), rho(:, p1:p2), 2);

            
            plot(t, Q)
            xlabel('Time [s]')
            ylabel('Charge [C cm-2]')
            xlim([t(1), t(end)])
        end

        function QVapp(sol, x1, x2)
            % Integrated charge density as a function of applied voltage
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            p1 = find(x<=x1);
            p1 = p1(end);
            p2 = find(x<=x2);
            p2 = p2(end);

            rho = dfana.calcrho(sol);
            Vapp = dfana.calcVapp(sol);
            Q = par.e*trapz(x(p1:p2), rho(:, p1:p2), 2);

            
            plot(Vapp, Q)
            xlabel('Vapp [V]')
            ylabel('Charge [C cm-2]')
            if Vapp(1) ~= Vapp(end)
                xlim([Vapp(1), Vapp(end)])
            end
        end

        function rhox(varargin)
            % Volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho(sol);

            
            dfplot.x2d(sol, x, {rho},{'\rho'},{'-'},'Charge density [cm-3]', tarr, xrange, 0, 0);
        end

        function rho_ionx(varargin)
            % Volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho_ion = dfana.calcrho_ion(sol);
            
            dfplot.x2d(sol, x, {rho_ion},{'\rho ion'},{'-'},'Charge density [cm-3]', tarr, xrange, 0, 0);
        end     
        
        function deltarhox(varargin)
            % The change in volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho(sol);
            deltarho = rho - rho(1,:);

            
            dfplot.x2d(sol, x, {deltarho},{'\Delta \rho'},{'-'},'Delta charge density [cm-3]', tarr, xrange, 0, 0);
        end

        function rhoxFxVx(varargin)
            % Three panel figure:
            % Volumetric charge density (rho), Field and potential as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            rho = dfana.calcrho(sol);
            F = dfana.calcF(sol);

            
            subplot(3, 1, 1)
            dfplot.x2d(sol, x, {rho},{'\rho'},{'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);

            subplot(3, 1, 2)
            dfplot.x2d(sol, x, {F},{'F'},{'-'},'Electric field [Vcm-1]', tarr, xrange, 0, 0);

            subplot(3, 1, 3)
            dfplot.x2d(sol, x, {V},{'V'},{'-'},'Electrostatic potential [V]', tarr, xrange, 0, 0);
        end

        function rhoxVx(varargin)
            % Three panel figure:
            % Volumetric charge density (rho), Field and potential as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            rho = dfana.calcrho(sol);

            
            subplot(2, 1, 1)
            dfplot.x2d(sol, x, {rho},{'\rho'},{'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);

            subplot(2, 1, 2)
            dfplot.x2d(sol, x, {-V},{'V'},{'-'},'-Electrostatic potential [V]', tarr, xrange, 0, 0);
        end

        function ELx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);

            
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'CB', 'VB'},...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        end

        function ELnpx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);

            
            subplot(2,1,1);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'CB', 'VB'},...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0);

            subplot(2,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, 'El carrier density [cm-3]', tarr, xrange, 0, 1);
        end

        function ELxnpxacx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);

            
            subplot(3,1,1);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'CB', 'VB'}, {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)

            subplot(3,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, 'El carrier density [cm-3]', tarr, xrange, 0, 1)

            
            subplot(3,1,3);
            dfplot.x2d(sol, x, {a, c}, {'a', 'c'}, {'-', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange, 0, 0)
        end

        function Vxacx(varargin)
            % Potential and ionic charges as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            
            subplot(2,1,1);
            dfplot.x2d(sol, x, {V}, {'V'},...
                {'-'}, 'Electro. potential [V]', tarr, xrange, 0, 0);

            subplot(2,1,2);
            dfplot.x2d(sol, x, {a, c}, {'a', 'c'},...
                {'-', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange , 0, 0);
        end

        function Vionxacx(varargin)
            % Electrostatic potential as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            Vion = dfana.calcVion(sol);
            Vel = V - Vion;

            
            subplot(2,1,1);
            dfplot.x2d(sol, x, {V, Vion, Vel}, {'V', 'Vion', 'Vel'},...
                {'--', '.', '-'}, 'Electro. potential [V]', tarr, xrange, 0, 0);

            subplot(2,1,2);
            dfplot.x2d(sol, x, {a, c}, {'a', 'c'},...
                {'-', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange , 0, 0);
        end

        function Fiont(sol, xpos)
            % Field contribution from ionic charge FION as a function of time at position XPOS
            [~,t,xmesh,~,~,~,~,~,~,~] = dfana.splitsol(sol);
            ppos = getpointpos(xpos, xmesh);
            Fion = dfana.calcFion(sol);
            
            
            plot(t, Fion(:,ppos))
            xlabel('Time')
            ylabel('Ion field [Vcm-1]')
        end

        function colourblocks(sol, yrange)
            par = sol.par;
            dcum0 = par.dcum0*1e7;   % Convert to nm

            for i =1:length(dcum0)-1
                v = [dcum0(i) yrange(2); dcum0(i+1) yrange(2); dcum0(i+1) yrange(1); dcum0(i) yrange(1)];   % vertices position
                f = [1 2 3 4];    % Faces
                if length(par.layer_colour) == length(dcum0)-1
                    j = i;
                else
                    j = i - ((length(par.layer_colour)-1)*floor(i/length(par.layer_colour)));
                end
                colour = par.layer_colour(j,:);
                patch('Faces',f,'Vertices',v,'FaceColor',colour, 'EdgeColor','none');%,'HandleVisibility','off')
            end

            hold on
        end

        function [sol, tarr, pointtype, xrange] = sortarg(args)

            if length(args) == 1
                sol = args{1};
                tarr = sol.t(size(sol.u,1));
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(args) == 2
                sol = args{1};
                tarr = args{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(args) == 3
                sol = args{1};
                tarr = args{2};
                xrange = args{3};
                pointtype = 't';
            end
        end

        function x2d(sol, xmesh, variables, legstr, linestyle, ylab, tarr, xrange, logx, logy)
            % SOL = solution structure
            % VARIABLES is an array containing the variables for plotting
            % LEGSTR is the legend string
            % YLAB = y-axis label
            % TARR- array of times
            % XRANGE - limits of the plot as a two element vector
            % LOGX, LOGY - switches for log axes
            ax = gca;
            if ishold(ax) == 0
                cla(ax);    % Clear current axis if held
            end
            par = sol.par;
            xnm = xmesh*1e7;

            vmin = min(min(cell2mat(variables)));
            vmax = max(max(cell2mat(variables)));

            if vmin == 0 && vmax == 0
                vmin = -1;
                vmax = 1;
            end
            vrange = vmax-vmin;
            if isempty(findobj(ax,'Type','patch'))
                switch logy
                    case 0
                        dfplot.colourblocks(sol, [vmin-(vrange*0.2), vmax+(vrange*0.2)]);
                    case 1
                        dfplot.colourblocks(sol, [0.1*vmin, 10*vmax]);
                end
            end

            vmin_tarr = zeros(length(tarr),length(variables));
            vmax_tarr = zeros(length(tarr),length(variables));
            h = zeros(1, length(variables));

            for i = 1:length(tarr)
                % find the time
                p1 = find(sol.t <= tarr(i));
                p1 = p1(end);
                for jj = 1:length(variables)
                    vtemp = variables{jj};

                    vmin_tarr(i,jj) = min(vtemp(p1, :));
                    vmax_tarr(i,jj) = max(vtemp(p1, :));

                    h(i,jj) = plot(xnm, variables{jj}(p1, :), char(linestyle(jj)));
                    hold on
                end
            end
            xlabel('Position [nm]')
            ylabel(ylab)
            if logy == 1
                set(gca, 'YScale','log');
            end
            if logx == 1
                set(gca, 'XScale','log');
            end
            if length(variables) == 1
                mystr = [];
                for i = 1:length(tarr)
                    mystr = [mystr, string(['t = ', num2str(tarr(i)), ' s'])];
                end
                lgd = legend(h, mystr);
            else
                lgd = legend(h(1,:), legstr);
            end
            lgd.FontSize = 12;
            xlim([xrange(1), xrange(2)])
            ymin = min(min(vmin_tarr));
            ymax = max(max(vmax_tarr));
            yrange = ymax - ymin;
            if ymin == 0 && ymax == 0
            else
                switch logy
                    case 0
                        if yrange == 0
                            ylim([ymin*0.9, ymax*1.1]);
                        else
                            ylim([ymin-(yrange*0.2), ymax+(yrange*0.2)]);
                        end
                    case 1
                        ylim([0.1*ymin, 10*ymax])
                end
            end
            set(gca, 'Layer', 'top')
            box on
            hold off
        end
    end
end
