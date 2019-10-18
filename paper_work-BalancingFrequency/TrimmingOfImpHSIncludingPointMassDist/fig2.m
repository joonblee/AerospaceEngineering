%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.2 Eigenfrequencies of the perfect shell                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%% variables %%%
syms phi;

%%% constants %%%
n=2; % mode number

E=210E9; rho=7800; mu=0.3;

aOh=[1:31]+9 % a/h
a=((15/(2*pi*rho))./(3./aOh+(1./aOh).^3/4)).^(1/3) % radius
%a=ones(1,31)
h=a./aOh % width

mass=4/3*pi*rho.*((a+h/2).^3-(a-h/2).^3)/2

%%% integrations %%%
U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,pi/2)
K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,pi/2)

%%% natural frequency of the perfect HS (omega0) %%%
omega0 = double(( n^2*(n^2-1)^2 * E*h.^2 ./ (3*(1+mu)*rho*a.^4) * U_phi/K_phi ) .^ .5) % angular freq
NatFreq0 = double(omega0 / (2*pi)) % linear freq

%%% figure %%%
figure(1);
grid on;
title('Eigenfrequencies of the perfect shell'); 
xlabel('a/h'); 
ylabel('Natural Frequency (=\omega_0/2\pi) [Hz]');

%%% plot (theory) %%%
figure(1);
hold on;
plot1(1)=plot(aOh,NatFreq0,'-','DisplayName','Analytical');

%%% plot (FEM) %%%
COMSOLx=[10;15;20;25;30;35;40]
COMSOLy=[1198;722;501;377;298;244;205]

figure(1);
hold on;
plot1(2)=plot(COMSOLx,COMSOLy,'x','DisplayName','FEM');

%%% save figure %%%
figure(1); legend(plot1(1:2)); legend boxoff;

saveas(figure(1),'figures/fig2.png');

%%% End of the code %%%
fprintf('# ----- END Calculation ----- # \n');
quit;
