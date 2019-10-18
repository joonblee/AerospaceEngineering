A = 8
I = 128/3/4
E = 10^7

q = 20

L1 = 30*sqrt(2)
L2 = 30

A1 = A %*sqrt(2)
I1 = I %*2*sqrt(2)

k1 = [A1*E/L1  0             0            -A1*E/L1 0             0;
      0        12*E*I1/L1^3  6*E*I1/L1^2  0        -12*E*I1/L1^3 6*E*I1/L1^2;
      0        6*E*I1/L1^2   4*E*I1/L1    0        -6*E*I1/L1^2  2*E*I1/L1;
	  -A1*E/L1 0             0            A1*E/L1  0             0;
	  0        -12*E*I1/L1^3 -6*E*I1/L1^2 0        12*E*I1/L1^3  -6*E*I1/L1^2;
	  0        6*E*I1/L1^2   2*E*I1/L1    0        -6*E*I1/L1^2  4*E*I1/L1]

R1 = 1/sqrt(2)*[1  1 0       0 0 0;
                -1 1 0       0 0 0;
				0  0 sqrt(2) 0 0 0;
				0  0 0       1 1 0;
				0  0 0       -1 1 0;
				0  0 0       0 0 sqrt(2)]

K1 = R1'*k1*R1


k2 = [A*E/L2  0            0           -A*E/L2 0            0;
      0       12*E*I/L2^3  6*E*I/L2^2  0       -12*E*I/L2^3 6*E*I/L2^2;
      0       6*E*I/L2^2   4*E*I/L2    0       -6*E*I/L2^2  2*E*I/L2;
      -A*E/L2 0            0           A*E/L2  0            0;
      0       -12*E*I/L2^3 -6*E*I/L2^2 0       12*E*I/L2^3  -6*E*I/L2^2;
      0       6*E*I/L2^2   2*E*I/L2    0       -6*E*I/L2^2  4*E*I/L2]

K2 = k2

Kstr = zeros(9);

for i = 1:6
  for j = 1:6
    Kstr(i,j) = Kstr(i,j) + K1(i,j);
  end
end

for i = 1:6
  for j = 1:6
    Kstr(i+3,j+3) = Kstr(i+3,j+3) + K2(i,j);
  end
end

Ffef = [0 0 0 0 -q*L2/2 q*L2^2/12 0 -q*L2/2 -q*L2^2/12]'
Fnode = [0 0 0 0 -1200 -1500 0 0 0]'
Fstr = Ffef+Fnode

Ksmall = Kstr(4:6,4:6)
Fsmall = Fstr(4:6)

U = inv(Ksmall)*Fsmall


f1 = K1 * 10^(-3) * [0 0 0 0.5137 -1.9978 0.0325]'
f2 = K2 * 10^(-3) * [0 0 0 0.5137 -1.9978 0.0325]' - [0 -300 -1500 0 -300 1500]'






quit();


