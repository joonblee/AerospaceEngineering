%%% Calculate connection factor(Lambda, L_{ijk}) & curvature(R_{ijkl}) in 4D

clear all;

syms t x y z c G M; 
q = [t x y z];

%{
g(:,:) = [c^2*(1+G*z/c^2)^2  0  0  0;
%g(:,:) = [              c^2  0  0  0;
                          0 -1  0  0;
		                  0  0 -1  0;
		                  0  0  0 -1]

%g_ = inv(g);

g_(:,:) = [1/c^2*1/(1+G*z/c^2)^2  0  0  0;
                               0 -1  0  0;
							   0  0 -1  0;
							   0  0  0 -1]
%}

% Schwarzschild metric (x:r, y:theta, z:phi)
g(:,:) = [1-2*G*M/(c^2*x)                    0    0               0;
                        0 -1/(1-2*G*M/(c^2*x))    0               0;
                        0  0                   -x^2               0;
                        0  0                      0 -x^2*(sin(y))^2]

g_(:,:) = [1/(1-2*G*M/(c^2*x))                0      0               0;
                             0 -1+2*G*M/(c^2*x)      0               0;
                             0                0 -1/x^2               0;
                             0                0      0 -1/(x^2*(sin(y))^2)]

% connection coefficient / christoffel symbol
L = sym('L',[length(g) length(g) length(g)]);

for i = 1:length(g)
  for j = 1:length(g)
    for k = 1:length(g)
      L(i,j,k) = 1/2*g_(i,1)*(diff(g(1,j),q(k))+diff(g(1,k),q(j))-diff(g(j,k),q(1)));
      for l = 2:length(g)
	    L(i,j,k) = L(i,j,k) + 1/2*g_(i,l)*(diff(g(l,j),q(k))+diff(g(l,k),q(j))-diff(g(j,k),q(l)));
      end
	end
  end
end

% curvature tensor
R = sym('R',[length(g) length(g) length(g) length(g)]);

for i = 1:length(g)
  for j = 1:length(g)
    for k = 1:length(g)
      for l = 1:length(g)
	    R(l,i,j,k) = diff(L(l,i,k),q(j))-diff(L(l,i,j),q(k))+L(1,i,k)*L(l,1,j)-L(1,i,j)*L(l,1,k);
		for m = 2:length(g)
		  R(l,i,j,k) = R(l,i,j,k)+L(m,i,k)*L(l,m,j)-L(m,i,j)*L(l,m,k);
		end
      end
	end
  end
end

disp(L);
disp(R);

