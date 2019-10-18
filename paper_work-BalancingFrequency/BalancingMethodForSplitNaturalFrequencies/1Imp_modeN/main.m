%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.3 Mode shape of the imperfect shell                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('--------------------- Start Code -----------------------\n');

syms phi;

E=210E9; rho=7800; mu=0.3;
M=5; phi0=0.5*pi;

hOa = 0.01; % a/h
a=((M*3/pi/rho)/(2+cos(phi0))/(1-cos(phi0))^2./(3*hOa+(hOa).^3/4)).^(1/3)
% 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg
h   = a*hOa; % width

n=[2:6]';

freq0_mass = [];
freq0_hole = [];
freqL_mass = [];
freqH_mass = [];
freqL_hole = [];
freqH_hole = [];


for pm = [1 -1]

m_i = pm*m_i

C_i = zeros(length(n),length(phi_i));
denominator = zeros(length(n),1);
numerator   = zeros(length(n),1);

for i = 1:length(theta_i)
  fprintf(['i =' i '\n']);
  C_i(:,i) = m_i(i)*tan(phi_i(i)/2).^(2*n).*(n+cos(phi_i(i))).^2
  numerator   = numerator + C_i(:,i).*sin(2*n*theta_i(i))
  denominator = denominator + C_i(:,i).*cos(2*n*theta_i(i))
end

%%% shift angle %%%
j=0;
zetaL = 1./(2*n) .* (atan(numerator./denominator)+j*pi) %( theta_i + pi ); % j = 0 = L
zetaH = 1./(2*n) .* (atan(numerator./denominator)+(j+1)*pi) %theta_i; % j = 1 = H

U_phi = zeros(length(n),1);
K_phi = zeros(length(n),1);

% integrations
for i=1:length(n)
  U_phi(i) = int( tan(phi/2)^(2*n(i)) / sin(phi)^3 ,phi,0,phi0);
  K_phi(i) = int( tan(phi/2)^(2*n(i)) * ((n(i)+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,phi0);
end

K_iL = zeros(length(n),1);
K_iH = zeros(length(n),1);
for i=1:length(phi_i)
  K_iL = K_iL + m_i(i) * tan(phi_i(i)/2).^(2*n) .* (sin(phi_i(i))^2 + (n+cos(phi_i(i))).^2 .* sin(n.*(theta_i(i)-zetaL)).^2 );
  K_iH = K_iH + m_i(i) * tan(phi_i(i)/2).^(2*n) .* (sin(phi_i(i))^2 + (n+cos(phi_i(i))).^2 .* sin(n.*(theta_i(i)-zetaH)).^2 );
end

% epsilon
epsilonK_L = K_iL./K_phi/(pi*rho*a^2*h)
epsilonK_H = K_iH./K_phi/(pi*rho*a^2*h)

% eigen frequency (angular freq, natural frequency)
omega0 = double(( n.^2.*(n.^2-1).^2 * E*h^2 ./ (3*(1+mu)*rho*a^4) .* U_phi./K_phi ).^.5)
omegaL = double(omega0./(1+epsilonK_L).^.5)
omegaH = double(omega0./(1+epsilonK_H).^.5)

% linear frequency
freq0 = double(omega0 / (2*pi))
freqL = double(omegaL / (2*pi))
freqH = double(omegaH / (2*pi))

if pm == 1
  freq0_mass = freq0;
  freqL_mass = freqL;
  freqH_mass = freqH;
else
  freqL_hole = freqL;
  freqH_hole = freqH;
end

end



fig=figure(1);

gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];

n
freq0_mass
freqL_mass
freqH_mass
freqL_hole
freqH_hole

p(1)=plot(n,freq0_mass,'-ko','DisplayName','perfect');
hold on;
p(2)=plot(n,freqL_mass,'-bx','DisplayName','1M (j = L)');
p(3)=plot(n,freqH_mass,'-+','color',[0.58 0 0.83],'DisplayName','1M (j = H)');
p(4)=plot(n,freqL_hole,'-rs','DisplayName','1H (j = L)');
p(5)=plot(n,freqH_hole,'-^','color',[1.0 0.5 0.0],'DisplayName','1H (j = H)');

set(gca,'XTick',2:6);
set(gca,'XTickLabel',{'2','3','4','5','6'});

xlabel('mode number');
ylabel('Frequency [Hz]');

legend(p(1:5),'Location','NorthWest'); legend boxoff;

% save plot
fsign  = '';
fphi   = '';
ftheta = '';

for i = 1:length(theta_i)
  fsign = [fsign SIGN_i(i)];
  fphi  = [fphi '-' num2str(rad2deg(phi_i(i)))];
  ftheta = [ftheta '-' num2str(rad2deg(theta_i(i)))];
end

saveas(fig,['figures/1imp_modeN.png']);
saveas(fig,['figures/1imp_modeN.pdf']);


fprintf('\n\n');

%%% clear figure %%%
clf(figure(1));

fprintf('# ----- END Calculation ----- # \n');
quit;
