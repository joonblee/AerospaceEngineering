%%% Calculate connection factor(Lambda, L_{ijk}) & curvature(R_{ijkl}) in 3D

clear all;

syms q L E I;

k(:,:) = E*I/L^3 * [ 12   6*L  -12    6*L;
                    6*L 4*L^2 -6*L  2*L^2;
					-12  -6*L   12   -6*L;
					6*L 2*L^2 -6*L -4*L^2]

inv(k)

f = q*L*[-3/20 -L/30 -7/20 L/20]'

U = inv(k)*f


