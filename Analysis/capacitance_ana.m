function [J_electronic, delta_t, dV_by_dT_across_points, C_as_function_V_across_points] = capacitance_ana(sol_CV)
par = sol_CV.par;
[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_CV);
J = dfana.calcJ(sol_CV, "sub");%Get current density from dfana
J_tot=J.tot;
J_disp=J.disp;
J_electronic=J.tot-J.disp;
%% Get delta t
delta_t=t(:,2)-t(:,1); %change in time interval 
%% Create loop to calculate change in potential
for i=1:length(t)-1
    for j=1:length(x)
       dV_by_dT_across_points(i,j)= (V(i,j)-V(i+1,j))/delta_t;  %this is change in potential at each place for different times
    end
end

%% Remove the first column since it does not contribute
dV_by_dT_across_points(:,1) = []; 
%% Remove final row since it does not contribute
J_disp(481,:)=[]; 
%% Calculate capacitance at each point
for i=1:length(t)-1
    for j=1:length(x)-1
        C_as_function_V_across_points(i,j) = J_disp(i,j)/dV_by_dT_across_points(i,j);%Get C=J(V)/(dV/dt)
    end
end
%% Plot capacitance as a function of position
x(:,1)=[]; %Remove the first point to since it doesnt come inside the calculation
figure(7464)
plot(x(100:320),C_as_function_V_across_points(6,(100:320))); 
end