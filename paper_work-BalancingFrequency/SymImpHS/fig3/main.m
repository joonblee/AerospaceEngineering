%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.4 Frequencies with respect to the ratio of imperfection                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%% constants %%%
n=2;

E=210E9; rho=7800; mu=0.3;

hOa = 0.01; % a/h
a   = ((15/(2*pi*rho))/(3*hOa+(hOa)^3/4))^(1/3); % radius
% 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg
h   = a*hOa; % width

M=5;

r= [0:.01:.1];

phi_i   = pi/2;
theta_i =    0;

%%% Run iteration %%%
% attached mass
m_i = r*M;

% variables
syms phi;

% shift angle
zetaL = 1/(2*n) * theta_i; % j = 0 = L
zetaH = 1/(2*n) * (theta_i + pi ); % j = 1 = H

% integrations
U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,pi/2);
K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,pi/2);

K_iL = tan(phi_i/2)^(2*n) * (sin(phi_i)^2 + (n+cos(phi_i))^2 * sin(n*(theta_i-zetaL))^2 );
K_iH = tan(phi_i/2)^(2*n) * (sin(phi_i)^2 + (n+cos(phi_i))^2 * sin(n*(theta_i-zetaH))^2 );

% epsilon
epsilonK_L = m_i*K_iL/K_phi/(pi*rho*a^2*h);
epsilonK_H = m_i*K_iH/K_phi/(pi*rho*a^2*h);

% eigen frequency (angular freq, natural frequency)
omega0 = double(( n^2*(n^2-1)^2 * E*h^2 ./ (3*(1+mu)*rho*a^4) * U_phi/K_phi )^.5);
omegaL = double(omega0./(1+epsilonK_L).^.5);
omegaH = double(omega0./(1+epsilonK_H).^.5);

% linear frequency
freq0 = double(omega0 / (2*pi)) * ones(1,length(r));
freqL = double(omegaL / (2*pi));
freqH = double(omegaH / (2*pi));

figure(1);
p(1)=plot(r,freq0,'color',[0 0.4470 0.7410],'LineStyle',':','LineWidth',2,'DisplayName','Perfect');
hold on;
p(2)=plot(r,freqL,'b','DisplayName','j = L mode');
p(3)=plot(r,freqH,'r--','DisplayName','j = H mode');

r
freq0
freqL
freqH

xlabel('m_{imp} / m_{HS}');
ylabel('Frequency [Hz]');

legend(p(1:3),'Location','best'); legend boxoff;

saveas(figure(1),'figures/fig3.png');

%%% End of the code %%%
fprintf('# ----- END Calculation ----- # \n');

quit;
