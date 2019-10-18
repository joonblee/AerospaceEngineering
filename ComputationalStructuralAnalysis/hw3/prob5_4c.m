%%% Calculate connection factor(Lambda, L_{ijk}) & curvature(R_{ijkl}) in 3D

clear all;

% functions = @root2d;
% c0 = [0,0];
% c = fsolve(functions,c0)

x=[0:0.01:1]';
% y_exact=-2/sin(1)*sin(x)+2*x;
y_fe1 = sqrt(168)*x.^2.*(1-x);
y_fe2 = 23.8094*(x.^2-x) - 13.8814*(x.^3-x);

figure(1);

grid on;
xlabel('x'); ylabel('y');

plot_y(1)=plot(x,y_fe1,'LineWidth',1,'DisplayName','y_{FE 1}');
hold on;
plot_y(2)=plot(x,y_fe2,'--','LineWidth',1,'DisplayName','y_{FE 2}');
%plot_y(3)=plot(x,y_fe2,':','LineWidth',2,'DisplayName','y_{FE 2}');

xlim([0 1]);

legend(plot_y(1:2),'Location','best'); legend boxoff;

saveas(figure(1),'fig5_4c.png');

quit;
