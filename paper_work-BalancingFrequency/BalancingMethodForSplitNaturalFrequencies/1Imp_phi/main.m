%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab independent run %
%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2; %% mode number

% imperfection %
MASS=0.1;

% ------------ %
phi0 = 2*pi/3;
THETA = 0;
PHI  = -phi0 + 2*[0:1800]/1800*phi0;

phi_i   = [PHI];
theta_i = [THETA];
m_i     = [MASS];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.3 Mode shape of the imperfect shell                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms phi;

E=210E9; rho=7800; mu=0.3;

M=5;

hOa = 0.01; % a/h
a=((M*3/pi/rho)/(2+cos(phi0))/(1-cos(phi0))^2./(3*hOa+(hOa).^3/4)).^(1/3);
%a   = ((15/(2*pi*rho))/(3*hOa+(hOa)^3/4))^(1/3); % radius for 5 kg HS
% 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg
h   = a*hOa; % width

Freq0=zeros(1,length(PHI));
Freq_Mass_L=zeros(1,length(PHI));
Freq_Mass_H=zeros(1,length(PHI));
Freq_Hole_L=zeros(1,length(PHI));
Freq_Hole_H=zeros(1,length(PHI));


for MassSign=[1 -1]

m_i = MassSign * m_i;

C_i         = zeros(length(theta_i),length(PHI));
denominator = zeros(1,length(PHI));
numerator   = zeros(1,length(PHI));

for i = 1:length(theta_i)
  % fprintf(['i =' i '\n']);
  C_i(i,:) = m_i(i)*tan(phi_i(i,:)/2).^(2*n).*(n+cos(phi_i(i,:))).^2;
  numerator   = numerator + C_i(i,:)*sin(2*n*theta_i(i));
  denominator = denominator + C_i(i,:)*cos(2*n*theta_i(i));
end

%%% shift angle %%%
j=0;
zetaL = 1/(2*n)*(atan(numerator./denominator)+j*pi); %( theta_i + pi ); % j = 0 = L
zetaH = 1/(2*n)*(atan(numerator./denominator)+(j+1)*pi); %theta_i; % j = 1 = H

% integrations
U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,phi0);
K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,phi0);


K_iL = zeros(1,length(PHI));
K_iH = zeros(1,length(PHI));

for i=1:length(theta_i)
  K_iL = K_iL + m_i(i) * tan(phi_i(i,:)/2).^(2*n) .* (sin(phi_i(i,:)).^2 + (n+cos(phi_i(i,:))).^2 .* sin(n*(theta_i(i)-zetaL)).^2 );
  K_iH = K_iH + m_i(i) * tan(phi_i(i,:)/2).^(2*n) .* (sin(phi_i(i,:)).^2 + (n+cos(phi_i(i,:))).^2 .* sin(n*(theta_i(i)-zetaH)).^2 );
end

%K_iL
%K_iH

% epsilon
epsilonK_L = K_iL/K_phi/(pi*rho*a^2*h);
epsilonK_H = K_iH/K_phi/(pi*rho*a^2*h);

% eigen frequency (angular freq, natural frequency)
omega0 = ones(1,length(PHI))*( n^2*(n^2-1)^2 * E*h^2 / (3*(1+mu)*rho*a^4) * U_phi/K_phi )^.5;
omegaL = omega0./(1+epsilonK_L).^.5;
omegaH = omega0./(1+epsilonK_H).^.5;

% linear frequency
freq0 = omega0 / (2*pi);
freqL = omegaL / (2*pi);
freqH = omegaH / (2*pi);

  if MassSign==1
    Freq0=freq0;
    Freq_Mass_L=freqL;
    Freq_Mass_H=freqH;
  else
    Freq_Hole_L=freqL;
    Freq_Hole_H=freqH;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%
% FEM (COMSOL) results %
%%%%%%%%%%%%%%%%%%%%%%%%
phi_FEM = [-2*pi/3, -pi/2, -pi/3, -pi/6, 0, pi/6, pi/3, pi/2, 2*pi/3];
Freq_Mass_perfect_FEM = 69.924;
Freq_Mass_L_FEM = [68.372, 69.686, 69.903, 69.922, 69.924, 69.922, 69.903, 69.686, 68.372];
Freq_Mass_H_FEM = [63.884, 68.693, 69.691, 69.870, 69.924, 69.870, 69.691, 68.693, 63.884];
Freq_Hole_perfect_FEM = 69.881;
Freq_Hole_L_FEM = [71.808, 70.192, 69.945, 69.892, 69.881, 69.892, 69.945, 70.192, 71.808];
Freq_Hole_H_FEM = [77.246, 71.230, 70.085, 69.923, 69.881, 69.923, 70.085, 71.230, 77.246];

%%%%%%%%%%
% Figure %
%%%%%%%%%%
fig=figure(1);

gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];

grid on;

fig;
hold on;
%p(0)=plot(PHI,Freq0,'k','LineWidth',1.1,'DisplayName','perfect'); %'color',[0 0.4470 0.7410],'Marker','o','DisplayName','Perfect');
p(1)=plot(PHI,Freq_Mass_L./Freq0,'b','LineStyle','--','LineWidth',.5,'DisplayName','1M, j = L (Analytic)');
p(2)=plot(phi_FEM,Freq_Mass_L_FEM/Freq_Mass_perfect_FEM,'bo','MarkerSize',3,'DisplayName','1M, j = L (FEM)');
p(3)=plot(PHI,Freq_Mass_H./Freq0,'color',[0.58 0 0.83],'LineStyle','--','LineWidth',.5,'DisplayName','1M, j = H (Analytic)');
p(4)=plot(phi_FEM,Freq_Mass_H_FEM/Freq_Mass_perfect_FEM,'d','MarkerSize',3,'color',[0.58 0 0.83],'DisplayName','1M, j = H (FEM)');
p(5)=plot(PHI,Freq_Hole_L./Freq0,'r','LineStyle',':','LineWidth',1.1,'DisplayName','1H, j = L (Analytic)');
p(6)=plot(phi_FEM,Freq_Hole_L_FEM/Freq_Mass_perfect_FEM,'rs','MarkerSize',3,'DisplayName','1H, j = L (FEM)');
p(7)=plot(PHI,Freq_Hole_H./Freq0,'color',[1.0 0.5 0.0],'LineStyle',':','LineWidth',1.1,'DisplayName','1H, j = H (Analytic)');
p(8)=plot(phi_FEM,Freq_Hole_H_FEM/Freq_Mass_perfect_FEM,'^','MarkerSize',3,'color',[1.0 0.5 0.0],'DisplayName','1H, j = H (FEM)');

set(gca,'XTick',-phi0:pi/6:phi0);
set(gca,'XTickLabel',{'-2\pi/3','-\pi/2','-\pi/3','-\pi/6','0','\pi/6','\pi/3','\pi/2','2\pi/3'});

%xlim([-1.571 1.571]);
xlim([-phi0-0.001 phi0+0.001]);

xlabel('\phi');
ylabel('Frequency Ratio ($r$)','Interpreter','Latex');

legend(p(1:8),'Location','NorthEast'); legend boxoff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End matlab independent run %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveas(fig,['1imp_phi.png']);
saveas(fig,['1imp_phi.pdf']);
%saveas(fig,['figures/fig9_n' num2str(n) '.png']);
%saveas(fig,['figures/fig9_n' num2str(n) '.pdf']);

fprintf('\n\n');

%%% clear figure %%%
clf(fig);

fprintf('# ----- END Calculation ----- # \n');
quit;

