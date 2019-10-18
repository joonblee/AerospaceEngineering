L = 24;
yMax = -0.3750;
theta1 = -0.0163;
v2 = -0.3750;
theta2 = -0.0147;

x=0

Sig = yMax*((6*x/L^2-4/L)*theta1 + (6/L^2-12*x/L^3)*v2 + (6*x/L^2-2/L)*theta2)

x=L

Sig = yMax*((6*x/L^2-4/L)*theta1 + (6/L^2-12*x/L^3)*v2 + (6*x/L^2-2/L)*theta2)

x=3/4*L

Sig = yMax*((6*x/L^2-4/L)*theta1 + (6/L^2-12*x/L^3)*v2 + (6*x/L^2-2/L)*theta2)

quit();


