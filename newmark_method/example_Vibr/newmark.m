function [t,d,v,a] = NewmarkMethod(beta,gamma,dt)

fprintf('\n# ----- Start Newmark Calculation ----- #\n\n');
% ----- step1: initial conditions ----- %
fprintf('----- initial conditions -----\n');

m = 10; e = 2; J = 20;
k1 = 50; k2 = 100;
%beta = 1/6; gamma = 1/2;
ti = 0.; tf = 10.0; %dt = 0.1;
t = ti:dt:tf;

M = [ m   -m*e     ;
     -m*e  m*e^2+J];

K = [ k1   0  ;
      0    k2];

P = [100 0]';

% d, v, a %
d = zeros(length(M),length(t)); v = zeros(length(M),length(t)); a = zeros(length(M),length(t));
for i = 1:length(t)
  a(:,i) = 999;
  if i == 1
    d(:,i) = 0; v(:,i) = 0;
  else
    d(:,i) = 999; v(:,i) = 999;
  end
end

% Ft %
Ft = zeros(length(M),length(t)); 
for i = 1: length(t)
  Ft(:,i) = P;
end

% print %
fprintf('Print initial variables\n');
M
K
P
fprintf('#####\n');
d
v
a
Ft

fprintf('------------------------------\n\n');

% --------------------------------------------------------------------------------------------- %
% ----- step2: initial acceleration ----- %
fprintf('----- initial acceleration -----\n');
a(:,1) = inv(M)*(Ft(:,1) - K*d(:,1));
fprintf('initial acceleration, a(1)');
a(:,1)
fprintf('--------------------------------\n\n');
% --------------------------------------- %

% --------------------------------------------------------------------------------------------- %
% ----- step3: run iteration ----- %
fprintf('----- run iteration -----\n');
K_ = M/(beta*dt^2)+K;
fprintf('K_ =');
K_

fprintf('_______________________________________________________________\n');
for i = 1:length(t)-1
  d(:,i+1) = inv(K_)*(Ft(:,i+1)+M/(beta*dt^2)*(d(:,i)+dt*v(:,i)+dt^2*(1/2-beta)*a(:,i)));
  a(:,i+1) = inv(M)*(Ft(:,i+1)-K*d(:,i+1));
  v(:,i+1) = v(:,i) + dt*((1-gamma)*a(:,i)+gamma*a(:,i+1));
end
fprintf('output\n');
d
v
a
fprintf('_______________________________________________________________\n\n');
% -------------------------------- %

fprintf('# ----- END Newmark Calculation ----- # \n');

