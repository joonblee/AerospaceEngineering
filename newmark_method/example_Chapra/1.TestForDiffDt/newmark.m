function [t,d,v,a] = NewmarkMethod(beta,gamma,dt)

fprintf('##### ============================================ #####\n');
fprintf('Start the code \n\n');

% ----- step1: initial conditions ----- %
fprintf('----- initial conditions -----\n');

m = 1.77; k = 70;
%beta = 1/6; gamma = 1/2;
ti = 0.; tf = 2.0; %dt = 0.1;
t = ti:dt:tf;

% d, v, a %
d = zeros(length(t),1); v = zeros(length(t),1); a = zeros(length(t),1);
for i = 1:length(t)
  a(i) = 999;
  if i == 1
    d(i) = 0; v(i) = 0;
  else
    d(i) = 999; v(i) = 999;
  end
end

% Ft %
Ft = zeros(length(t),1); 
for i = 1: length(t)
  if t(i) < 0.25
    Ft(i) = 100 - 200 * t(i);
  else
    Ft(i) = 50 - 50/1.75 * (t(i)-0.25);
  end
end

% print %
fprintf('m = %.2f, k = %d                \n',m,k);
fprintf('beta = %s, gamma = %s           \n',strtrim(rats(beta)),strtrim(rats(gamma)));
fprintf('ti = %.2f, tf = %.2f, dt = %.2f \n',ti,tf,dt);
fprintf('d0 = %f, v0 = %f                \n',d(1),v(1));
fprintf('------------------------------\n\n\n');

% --------------------------------------------------------------------------------------------- %
% ----- step2: initial acceleration ----- %
fprintf('----- initial acceleration -----\n');
a(1) = (Ft(1) - k*d(1))/m;
fprintf('a(t=0) = %.2f       \n',a(1));
fprintf('--------------------------------\n\n\n');
% --------------------------------------- %

% --------------------------------------------------------------------------------------------- %
% ----- step3: run iteration ----- %
fprintf('----- run iteration -----\n');
K = m/(beta*dt^2)+k;
fprintf('k'' = %s [m/(beta*dt^2)+k]\n',K);

fprintf('_______________________________________________________________\n');
fprintf('t = 0    | Ft = %7.2f | d = %6.2f | v = %6.2f | a = %7.2f\n',Ft(1),d(1),v(1),a(1));
for i = 1:length(t)-1
  d(i+1) = (Ft(i+1)+m/(beta*dt^2)*(d(i)+dt*v(i)+dt^2*(1/2-beta)*a(i)))/K;
  a(i+1) = (Ft(i+1)-k*d(i+1))/m;
  v(i+1) = v(i) + dt*((1-gamma)*a(i)+gamma*a(i+1));
  fprintf('    %.2f |      %7.2f |     %6.2f |     %6.2f |     %7.2f\n',t(i+1),Ft(i+1),d(i+1),v(i+1),a(i+1));
end
fprintf('_______________________________________________________________\n');
fprintf('-------------------------\n\n\n');
% -------------------------------- %

fprintf('\n ##### END Newmark Calculation ##### \n\n');

