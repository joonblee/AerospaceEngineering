%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.3 Mode shape of the imperfect shell                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('--------------------- main.m -----------------------\n');

syms phi;

E=210E9; rho=7800; mu=0.3;

hOa = 0.01; % a/h
a   = ((15/(2*pi*rho))/(3*hOa+(hOa)^3/4))^(1/3); % radius for 5 kg HS
% 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg
h   = a*hOa; % width


C_i         = zeros(length(theta_i),length(PHI1));
denominator = zeros(1,length(PHI1));
numerator   = zeros(1,length(PHI1));

for i = 1:length(theta_i)
  % fprintf(['i =' i '\n']);
  C_i(i,:) = m_i(i)*tan(phi_i(i,:)/2).^(2*n).*(n+cos(phi_i(i,:))).^2
  numerator   = numerator + C_i(i,:)*sin(2*n*theta_i(i))
  denominator = denominator + C_i(i,:)*cos(2*n*theta_i(i))
end

%%% shift angle %%%
j=0;
zetaL = 1/(2*n)*(atan(numerator./denominator)+j*pi) %( theta_i + pi ); % j = 0 = L
zetaH = 1/(2*n)*(atan(numerator./denominator)+(j+1)*pi) %theta_i; % j = 1 = H

% integrations
U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,pi/2);
K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,pi/2);


K_iL = zeros(1,length(PHI1));
K_iH = zeros(1,length(PHI1));

for i=1:length(theta_i)
  K_iL = K_iL + m_i(i) * tan(phi_i(i,:)/2).^(2*n) .* (sin(phi_i(i,:)).^2 + (n+cos(phi_i(i,:))).^2 .* sin(n*(theta_i(i)-zetaL)).^2 )
  K_iH = K_iH + m_i(i) * tan(phi_i(i,:)/2).^(2*n) .* (sin(phi_i(i,:)).^2 + (n+cos(phi_i(i,:))).^2 .* sin(n*(theta_i(i)-zetaH)).^2 )
end

%K_iL
%K_iH

% epsilon
epsilonK_L = K_iL/K_phi/(pi*rho*a^2*h)
epsilonK_H = K_iH/K_phi/(pi*rho*a^2*h)

% eigen frequency (angular freq, natural frequency)
omega0 = ones(1,length(PHI1))*( n^2*(n^2-1)^2 * E*h^2 / (3*(1+mu)*rho*a^4) * U_phi/K_phi )^.5
omegaL = omega0./(1+epsilonK_L).^.5
omegaH = omega0./(1+epsilonK_H).^.5

% linear frequency
freq0 = omega0 / (2*pi)
freqL = omegaL / (2*pi)
freqH = omegaH / (2*pi)

figure(1);
grid on;

figure(1);
hold on;
p(1)=plot(PHI2,freq0,'color',[0 0.4470 0.7410],'DisplayName','Perfect'); %'color',[0 0.4470 0.7410],'Marker','o','DisplayName','Perfect');
p(2)=plot(PHI2,freqL,'b','DisplayName','j = L mode');
p(3)=plot(PHI2,freqH,'r--','DisplayName','j = H mode');

set(gca,'XTick',-pi/2:pi/2:pi/2);
set(gca,'XTickLabel',{'-\pi/2','0','\pi/2'});

xlim([-1.571 1.571]);

xlabel('\Delta\phi');
ylabel('Frequency [Hz]');

legend(p(1:3),'Location','NorthEast'); legend boxoff;

fsign = '';
for i = 1:length(SIGN_i)
  fsign = [fsign SIGN_i(i)];
end

saveas(figure(1),['figures/fig9_n' num2str(n) fsign '.png']);


fprintf('\n\n');

%%% clear figure %%%
clf(figure(1));

fprintf('# ----- END Calculation ----- # \n');
quit;

