clc;
fprintf('\n### ----- Start run.m ----- ###\n');

figure(1);
line([0 2],[0 0],'color','black');
grid on;
title('t vs d'); xlabel('time (s)'); ylabel('distance (m)');

figure(2);
line([0 2],[0 0],'color','black');
grid on;
title('t vs v'); xlabel('time (s)'); ylabel('velocity (m/s)');

figure(3);
line([0 2],[0 0],'color','black');
grid on;
title('t vs a'); xlabel('time (s)'); ylabel('acceleration (m/s^2)');

plot_style = ["b","rs","m^","kx","g+","c*","y."];

loop = 1; Nth_loop = 1;

while loop
  fprintf('\n# ----- Start %d th loop ----- # \n\n',Nth_loop);
  clear beta; clear gamma; clear dt;
  clear t; clear d; clear v; clear a;

  str = 'Input beta: ';
  beta=input(str);
  str = 'Input gamma: ';
  gamma=input(str);
  str = 'Input dt: ';
  dt=input(str);

  [t,d,v,a] = newmark(beta,gamma,dt);

  % ----- plotting ----- %
  figure(1);
  hold on;
  plot_d(Nth_loop)=scatter(t,d,plot_style(Nth_loop),'DisplayName',['\beta=',num2str(beta),', \gamma=',num2str(gamma),', \Deltat=',num2str(dt)]);

  figure(2);
  hold on;
  plot_v(Nth_loop)=scatter(t,v,plot_style(Nth_loop),'DisplayName',['\beta=',num2str(beta),', \gamma=',num2str(gamma),', \Deltat=',num2str(dt)]);

  figure(3);
  hold on;
  plot_a(Nth_loop)=scatter(t,a,plot_style(Nth_loop),'DisplayName',['\beta=',num2str(beta),', \gamma=',num2str(gamma),', \Deltat=',num2str(dt)]);

  Nth_loop = Nth_loop + 1;

  str = '\nDo you want to repeat a loop? [yes:1/no:0] ';
  loop = input(str);
  fprintf('\n');
end;

figure(1); legend(plot_d(1:Nth_loop-1)); legend boxoff;
figure(2); legend(plot_v(1:Nth_loop-1)); legend boxoff;
figure(3); legend(plot_a(1:Nth_loop-1)); legend boxoff;

saveas(figure(1),'distance.png');
saveas(figure(2),'velocity.png');
saveas(figure(3),'acceleration.png');

fprintf('# ===== SUCCESS ===== #\n\n');
quit;
