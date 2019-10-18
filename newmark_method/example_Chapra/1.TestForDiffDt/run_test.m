
beta = 1/6; gamma = 1/2;
dt = 0.1;
[t1,d1,v1,a1] = newmark(beta,gamma,dt);

dt = 0.05;
[t2,d2,v2,a2] = newmark(beta,gamma,dt);

dt = 0.01;
[t3,d3,v3,a3] = newmark(beta,gamma,dt);

% ----- plotting ----- %
figure(1);

line([0 2],[0 0],'color','black');
grid on;
hold on;

Plot_d1 = scatter(t1,d1,'b');
Plot_d2 = scatter(t2,d2,'rs');
Plot_d3 = scatter(t3,d3,'m.');
%axis(xmin xmax ymin ymax);
%plot(t,d);

title('t vs d');
xlabel('time (s)');
ylabel('distance (m)');
legend([Plot_d1 Plot_d2 Plot_d3],{'dt=0.1','dt=0.05','dt=0.01'});
hold off;

saveas(figure(1),'distance.png');

% ----- %
figure(2);
line([0 2],[0 0],'color','black');
grid on;
hold on;

Plot_v1 = scatter(t1,v1,'b');
Plot_v2 = scatter(t2,v2,'rs');
Plot_v3 = scatter(t3,v3,'m.');

title('t vs v');
xlabel('time (s)');
ylabel('velocity (m/s)');
legend([Plot_v1 Plot_v2 Plot_v3],{'dt=0.1','dt=0.05','dt=0.01'});
hold off;

saveas(figure(2),'velocity.png');

% ----- %
figure(3);

line([0 2],[0 0],'color','black');
grid on;
hold on;

Plot_a1 = scatter(t1,a1,'b');
Plot_a2 = scatter(t2,a2,'rs');
Plot_a3 = scatter(t3,a3,'m.');

title('t vs a');
xlabel('time (s)');
ylabel('accel (m/s^2)');
legend([Plot_a1 Plot_a2 Plot_a3],{'dt=0.1','dt=0.05','dt=0.01'});
hold off;

saveas(figure(3),'acceleration.png');
% ------------------------------------------ %

fprintf('# ===== SUCCESS ===== #\n\n');
quit;
