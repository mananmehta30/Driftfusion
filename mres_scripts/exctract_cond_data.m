% script to get data points from conductance plot

fig = openfig('conductances0V_Ag-MAPI-Al.fig');
%%

A = 1e-4; % 0.01mm-2 = 1e-4 cm-2

fig = gcf;

dataObjs = findobj(fig,'-property','YData');

x_both = dataObjs(2).XData;
y_both = dataObjs(2).YData;
R_both = 1./y_both/A;

x_l = dataObjs(3).XData;
y_l = dataObjs(3).YData;
R_l = 1./y_l/A;

x_r = dataObjs(4).XData;
y_r = dataObjs(4).YData;
R_r = 1./y_r/A;

x_no = dataObjs(5).XData;
y_no = dataObjs(5).YData;
R_no = 1./y_no/A;

x_el = dataObjs(6).XData;
y_el = dataObjs(6).YData;
R_el = 1./y_el/A;

R


figure
plot(x_el,R_el)
hold on
plot(x_no,R_no)
plot(x_r,R_r)
plot(x_l,R_l)
plot(x_both,R_both)
hold off
set(gca,'yscale','log')