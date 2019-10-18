%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.3 Mode shape of the imperfect shell                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('--------------------- main.m -----------------------\n');

syms phi;

n=2;
E=210E9; rho=7800; mu=0.3;

m_i     = [0.1; -0.05; 0.07]
phi_i   = [pi/4;  pi/2; 2*pi/3]
theta_i = [0;     pi/2; 4*pi/3]

M=5
phi0=2*pi/3

hOa = 0.01; % a/h
a=((M*3/pi/rho)/(2+cos(phi0))/(1-cos(phi0))^2./(3*hOa+(hOa).^3/4)).^(1/3)
% 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg
h   = a*hOa; % width


C_i         = 0;
denominator = 0;
numerator   = 0;

for i = 1:length(theta_i)
  % fprintf(['i =' i '\n']);
  C_i = m_i(i)*tan(phi_i(i)/2)^(2*n)*(n+cos(phi_i(i)))^2;
  numerator   = numerator + C_i*sin(2*n*theta_i(i));
  denominator = denominator + C_i*cos(2*n*theta_i(i));
end

%%% shift angle %%%
j=0;
zetaL = 1/(2*n)*(atan(numerator/denominator)+j*pi) %( theta_i + pi ); % j = 0 = L
zetaH = 1/(2*n)*(atan(numerator/denominator)+(j+1)*pi) %theta_i; % j = 1 = H

% integrations
U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,phi0);
K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,phi0);


K_L = 0;
K_H = 0;

for i=1:length(theta_i)
  K_L = K_L + m_i(i) * tan(phi_i(i)/2)^(2*n) * (sin(phi_i(i))^2 + (n+cos(phi_i(i)))^2 * sin(n*(theta_i(i)-zetaL))^2 )
  K_H = K_H + m_i(i) * tan(phi_i(i)/2)^(2*n) * (sin(phi_i(i))^2 + (n+cos(phi_i(i)))^2 * sin(n*(theta_i(i)-zetaH))^2 )
end

% epsilon
epsilonK_L = K_L/K_phi/(pi*rho*a^2*h)
epsilonK_H = K_H/K_phi/(pi*rho*a^2*h)

% eigen frequency (angular freq, natural frequency)
omega0 = ( n^2*(n^2-1)^2 * E*h^2 / (3*(1+mu)*rho*a^4) * U_phi/K_phi )^.5
omegaL = omega0/(1+epsilonK_L)^.5
omegaH = omega0/(1+epsilonK_H)^.5


phi_b = [0:200]/200*phi0;
m_b = a^2*h*pi*rho*n^2*K_phi*(omegaH^2-omegaL^2)/((1+n^2)*omegaH^2-omegaL^2)./(tan(phi_b/2).^(2*n).*(n+cos(phi_b)).^2);

kk=[0:2*n-1];
theta_b = zetaL + kk*pi/n;

PHIb   = phi_b' * ones(1,length(kk));
Mb     = m_b'   * ones(1,length(kk));
THETAb = ones(1,length(phi_b))' * theta_b;

PHIb = PHIb(:);
Mb   = Mb(:);
THETAb = THETAb(:);

eraser = zeros(1,length(Mb));
for i=[1:length(Mb)]
  if Mb(i) > M*0.2
    eraser(i)=1;
  elseif -Mb(i) > M*0.2
    eraser(i)=1;
  end
end

for i=[length(Mb):-1:1]
  if eraser(i)
    PHIb(i)=[];
	Mb(i)=[];
	THETAb(i)=[];
  end
end


PHIb = [PHIb; PHIb];
Mb   = [Mb;   -Mb];
THETAb = [THETAb; THETAb+pi/2/n];


fig=figure(1);

gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];

p(1)=polarscatter(THETAb,PHIb,5,Mb,'filled'); %'color',[0 0.4470 0.7410],'Marker','o','DisplayName','Perfect');
colorbar


pax = gca;
pax.RMinorGrid = 'on';

pax.ThetaTickLabel = [{'\theta = 0^{\circ}','\theta = 30^{\circ}','\theta = 60^{\circ}','\theta = 90^{\circ}','\theta = 120^{\circ}','\theta = 150^{\circ}','\theta = 180^{\circ}','\theta = 210^{\circ}','\theta = 240^{\circ}','\theta = 270^{\circ}','\theta = 300^{\circ}','\theta = 330^{\circ}'}] ;

rlim([0 phi0]);
rticks([0 pi/6 pi/3 pi/2 2*pi/3]);
%rticklabels({'\phi = 0','\phi = \pi/6','\phi = \pi/3','\phi = \pi/2','\phi = 2\pi/3'});
rticklabels({'\phi = 0^{\circ}','\phi = 30^{\circ}','\phi = 60^{\circ}','\phi = 90^{\circ}','\phi = 120^{\circ}'});


saveas(fig,'figures/fig.png');
saveas(fig,'figures/fig.pdf');


fprintf('\n\n');

%%% clear figure %%%
clf(figure(1));

fprintf('# ----- END Calculation ----- # \n');
quit;

