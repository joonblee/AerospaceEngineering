%%% Calculate connection factor(Lambda, L_{ijk}) & curvature(R_{ijkl}) in 3D

clear all;

syms r theta phi;
x = [theta phi r];

g(:,:) = [r^2 0 0;
          0 r^2*sin(theta)^2 0;
		  0 0 1]
% g(:,:) = [1 0 0; 0 1 0; 0 0 1];

g_ = inv(g)

% connection coefficient / christoffel symbol
L = sym('L',[length(g) length(g) length(g)]);

for i = 1:length(g)
  for j = 1:length(g)
    for k = 1:length(g)
      L(i,j,k) = 1/2*g_(i,1)*(diff(g(1,j),x(k))+diff(g(1,k),x(j))-diff(g(j,k),x(1)));
      for l = 2:length(g)
	    L(i,j,k) = L(i,j,k) + 1/2*g_(i,l)*(diff(g(l,j),x(k))+diff(g(l,k),x(j))-diff(g(j,k),x(l)));
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
	    R(l,i,j,k) = diff(L(l,i,k),x(j))-diff(L(l,i,j),x(k))+L(1,i,k)*L(l,1,j)-L(1,i,j)*L(l,1,k);
		for m = 2:length(g)
		  R(l,i,j,k) = R(l,i,j,k)+L(m,i,k)*L(l,m,j)-L(m,i,j)*L(l,m,k);
		end
      end
	end
  end
end

disp(L);
disp(R);


