%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.4 Frequencies with respect to the ratio of imperfection                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%% constants %%%
n=2;

E=210E9; rho=7800; mu=0.3;

aOh = 40; % a/h
a   = ((15/(2*pi*rho))/(3/aOh+(1/aOh)^3/4))^(1/3); % radius
h   = a/aOh; % width


M = 4/3*pi*rho.*((a+h/2).^3-(a-h/2).^3)/2;

r= [0:.001:.1];

phi_i   = pi/2;
theta_i =    0;

%%% Run iteration %%%
  % attached mass
  m_i = r*M;

  % variables
  syms phi;

  % shift angle
  zetaL = 1/(2*n) * ( theta_i + pi ); % j = 1 = L
  zetaH = 1/(2*n) * theta_i; % j = 2 = H

  % integrations
  U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,pi/2);
  K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,pi/2);

  K_iL = tan(phi_i/2)^(2*n) * (sin(phi_i)^2 + (n+cos(phi_i))^2 * sin(n*(theta_i-zetaL))^2 );
  K_iH = tan(phi_i/2)^(2*n) * (sin(phi_i)^2 + (n+cos(phi_i))^2 * sin(n*(theta_i-zetaH))^2 );

  % epsilon
  epsilonK_L = m_i*K_iL/K_phi/(pi*rho*a^2*h);
  epsilonK_H = m_i*K_iH/K_phi/(pi*rho*a^2*h);

  % eigen frequency (angular freq)
  omega0 = double(( n^2*(n^2-1)^2 * E*h^2 ./ (3*(1+mu)*rho*a^4) * U_phi/K_phi )^.5);
  omegaL = double(omega0./(1+epsilonK_L).^.5);
  omegaH = double(omega0./(1+epsilonK_H).^.5);

  % linear frequency
  NatFreq0 = double(omega0 / (2*pi)) * ones(1,length(r));
  NatFreqL = double(omegaL / (2*pi));
  NatFreqH = double(omegaH / (2*pi));
%end;

figure(1);
p(1)=plot(r,NatFreq0,'b','DisplayName','Perfect');
hold on;
p(2)=plot(r,NatFreqH,'r--','DisplayName','j = H mode');
p(3)=plot(r,NatFreqL,'m:','DisplayName','j = L mode');

xlabel('m (attached) / M (shell)');
ylabel('Natural frequency [Hz]');

legend(p(1:3),'Location','best'); legend boxoff;

saveas(figure(1),'figures/fig4.png');

%%% End of the code %%%
fprintf('# ----- END Calculation ----- # \n');

quit;
