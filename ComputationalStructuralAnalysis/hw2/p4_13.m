%%% Calculate connection factor(Lambda, L_{ijk}) & curvature(R_{ijkl}) in 3D

clear all;

k1(:,:) = 482.25/4 * [2304 -144 1152;
                   -144 12 -144;
				   1152 -144 2304]
k2(:,:) = [0 0 0;
           0 400 0;
		   0 0 0]

K = k1 + k2


f = [-225 -168.75 675]'

U = inv(K)*f

f1 = k1 * U
f2 = k2 * U

quit()

