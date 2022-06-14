figure
plot(Vsucc, JVgrad0_el_only, '.k', 'MarkerSize', 10)
hold on
plot(Vsucc, JVgrad0_no_sc,'+b', 'MarkerSize', 10)
plot(Vsucc, JVgrad0_sc1,'xr', 'MarkerSize', 10)
hold off

title('Conductance around 0V for different ionic behaviour')
legend('el only','no sc','sc_r = 1')
set(gca,'yscale','log')
xlabel('prebias voltage [V]')
ylabel('conductance [\Omega^{-1}cm^{-2}]')
%%
% Script to take excerpts and variables from B7 scripts and run them
% separately

Vi = 1;
Vp = -1;

ai = 1;

% sol_pre(Vi) = doJV(sol_eq(ai).ion, 1e-10, 101, 0, 1, 0, Vp, 1);
sol_pre(Vi) = jumptoV(sol_eq.ion, Vp, 500, 1, 0, 0, 0);

% sol_0V(Vi) = doJV(sol_pre(Vi).dk.f, 0.1, 101, 0, 0, Vp, 0, 1);
sol_0V(Vi) = jumptoV(sol_pre(Vi), 0, 500, 0, 0, 1, 0);


sol_CV_fr(Vi) = doCV(sol_0V(Vi), 0, 0,5, -5, 0.01, 1, 4001);
%% erxtract data from plot

fig = gcf;
fullsize = [1:points];
dataObjs = flipud(findobj(fig,'-property','YData'));

X = [];
Y = [];
for ii = 1:length(dataObjs)
    x = dataObjs(ii).XData;
%     x_full = [x , nan(1,numel(fullsize)-numel(x))];
%     X(ii,:) = x_full;
    
    y = dataObjs(ii).YData;
%     y_full = [y , nan(1,numel(fullsize)-numel(y))];
%     Y(ii,:) = y_full;

    JV_gradient = gradient(y,x);
    JV_gradient0(ii) = JV_gradient(1);
end


Voltages = [-3.7,-3.3,-2.3,-2.2,-2,-1.1,-.5,-.3,-.1,0,0.1,0.2,0.3,0.4,...
    0.5,0.6,0.8,0.9,1,1.1,1.2,1.5,1.6,1.7,1.8,2,2.1,2.2,2.3,2.4,2.5,2.7,...
    3,3.2,3.3,3.4,3.5,3.6,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5];

%%
figure
plot(Voltages(JV_gradient0 > 0),JV_gradient0(JV_gradient0 > 0), 'xb')
hold on
plot(Voltages(JV_gradient0 < 0),-JV_gradient0(JV_gradient0 < 0), 'xr')
hold off
set(gca, 'yscale', 'log')
legend('positive slope','negative slope')
title('Resistance as function of prebias, anions false')
xlabel('prebias voltage [V]')
ylabel('resistance [\Omega]')
