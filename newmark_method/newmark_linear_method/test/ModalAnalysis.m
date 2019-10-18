function [M_,K_,C_,P_,phi]= ModalAnalysis(M,K,C,P,lambda)
% Modal Analysis of the System 
%-------------------------------------------------------------------------
% Code written by :Siva Srinivas Kolukula                                 |
%                  Senior Research Fellow                                 |
%                  Structural Mechanics Laboratory                        |
%                  Indira Gandhi Center for Atomic Research               |
%                  INDIA                                                  |
% E-mail : allwayzitzme@gmail.com                                         |
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Purpose : To do the Modal Analysis of the system
%
% Synopsis : [M,K,C,P] = ModalAnalysis(M,K,C,P,lambda]
%
% Variable Description :
%           M - Mass Matrix
%           K - Stiffness Matrix
%           C - Damping Matrix
%           P - Force MAtrix
%           lambda - number of modes to be considered
%           phi - Reduced Normal modes of the system
%--------------------------------------------------------------------------

fprintf('\n\n__________________________________________\n');
fprintf(' ### ----- Modal Analysis ----- ### \n\n');
 % Modal Analysis 
fprintf('\n [V,D]=eig(K,M)\n');
fprintf('V: right eigenvectors, D: diagonal matrix of eigenvalues (KV = DMV)');
[V,D]=eig(K,M)

fprintf('\n%%%%%\n');
fprintf(' [W,k]=sort(diag(D))\n');
fprintf('_k_ is the same size as _diag(D)_ and describes the arrangement of the elements of _diag(D)_ into _W_ along the sorted dimension.\n');
fprintf('diag(D) =');
diag(D)
[W,k]=sort(diag(D))

fprintf('\n%%%%%\n');
fprintf(' V=V(:,k) \n');
fprintf('k_th column of matrix V.');
V=V(:,k)

fprintf('\n%%%%%\n');
fprintf(' Factor=diag(V''*M*V)\n');
fprintf('V''MV =');
V'*M*V
Factor=diag(V'*M*V)

fprintf('\n%%%%%\n');
fprintf(' Phi=V*inv(sqrt(diag(Factor)))\n');
fprintf('diag(Factor) =');
diag(Factor)
fprintf('sqrt(diag(Factor)) =');
sqrt(diag(Factor))
fprintf('inv(sqrt(diag(Factor))) =');
inv(sqrt(diag(Factor)))
Phi=V*inv(sqrt(diag(Factor)))

fprintf('\n%%%%%\n');
% selecting only first two Natural Frequencies and Mode shapes
lambda

fprintf('phi = Phi(:,1:lambda)');
phi = Phi(:,1:lambda)

fprintf('\n%%%%%\n');
fprintf('M_ = phi''Mphi, C_ = phi''Cphi, K_ = phi''Kphi, P_ = phi''p');
% Reducing the Degree's of Freedom of M, K and P
M_ = phi'*M*phi
C_ = phi'*C*phi
K_ = phi'*K*phi
P_ = phi'*P

fprintf('__________________________________________\n\n\n');

fprintf('### ----- Natural Frequencies, w ----- ### \n');

fprintf('\n Omega=diag(sqrt(Phi''*K*Phi))\n');
fprintf('Phi''KPhi =');
Phi'*K*Phi
fprintf('sqrt(Phi''KPhi) =');
sqrt(Phi'*K*Phi)
Omega=diag(sqrt(Phi'*K*Phi))

fprintf('\n w = diag(Omega(1:lambda))');
w = diag(Omega(1:lambda))

fprintf('__________________________________________\n\n\n');

quit;


