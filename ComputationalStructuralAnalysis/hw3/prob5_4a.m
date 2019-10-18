%%% Calculate connection factor(Lambda, L_{ijk}) & curvature(R_{ijkl}) in 3D

clear all;

x=[0:0.01:1]';
y_exact=-2/sin(1)*sin(x)+2*x;
y_fe1=-0.5*x.*(1-x.^2);
y_fe2=0.04336*(x.^2-x)+0.3415*(x.^3-x);

figure(1);

grid on;
xlabel('x'); ylabel('y');

plot_y(1)=plot(x,y_exact,'LineWidth',1,'DisplayName','y_{exact}');
hold on;
plot_y(2)=plot(x,y_fe1,'--','LineWidth',1,'DisplayName','y_{FE 1}');
plot_y(3)=plot(x,y_fe2,':','LineWidth',2,'DisplayName','y_{FE 2}');

xlim([0 1]);

legend(plot_y(1:3),'Location','best'); legend boxoff;

saveas(figure(1),'fig5_4a.png');

quit;
