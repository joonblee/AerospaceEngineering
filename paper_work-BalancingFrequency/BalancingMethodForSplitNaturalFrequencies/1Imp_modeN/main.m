%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab independent run %
%%%%%%%%%%%%%%%%%%%%%%%%%%  
m_i     = [0.1];
phi_i   = [pi/2];
theta_i = [0];

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
a=((M*3/pi/rho)/(2+cos(phi0))/(1-cos(phi0))^2./(3*hOa+(hOa).^3/4)).^(1/3);
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
  m_i = pm*m_i;

  C_i = zeros(length(n),length(phi_i));
  denominator = zeros(length(n),1);
  numerator   = zeros(length(n),1);

  for i = 1:length(theta_i)
    C_i(:,i) = m_i(i)*tan(phi_i(i)/2).^(2*n).*(n+cos(phi_i(i))).^2;
    numerator   = numerator + C_i(:,i).*sin(2*n*theta_i(i));
    denominator = denominator + C_i(:,i).*cos(2*n*theta_i(i));
  end

  %%% shift angle %%%
  j=0;
  zetaL = 1./(2*n) .* (atan(numerator./denominator)+j*pi); %( theta_i + pi ); % j = 0 = L
  zetaH = 1./(2*n) .* (atan(numerator./denominator)+(j+1)*pi); %theta_i; % j = 1 = H

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
  epsilonK_L = K_iL./K_phi/(pi*rho*a^2*h);
  epsilonK_H = K_iH./K_phi/(pi*rho*a^2*h);

  % eigen frequency (angular freq, natural frequency)
  omega0 = double(( n.^2.*(n.^2-1).^2 * E*h^2 ./ (3*(1+mu)*rho*a^4) .* U_phi./K_phi ).^.5);
  omegaL = double(omega0./(1+epsilonK_L).^.5);
  omegaH = double(omega0./(1+epsilonK_H).^.5);

  % linear frequency
  freq0 = double(omega0 / (2*pi));
  freqL = double(omegaL / (2*pi));
  freqH = double(omegaH / (2*pi));

  if pm == 1
    freq0_mass = freq0;
    freqL_mass = freqL;
    freqH_mass = freqH;
  else
    freqL_hole = freqL;
    freqH_hole = freqH;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%
% FEM (COMSOL) results %
%%%%%%%%%%%%%%%%%%%%%%%%
freq0_FEM      = [49.06 , 134.82, 254.60, 404.30, 580.90];
freqL_mass_FEM = [48.462, 133.51, 252.59, 401.72, 577.70];
freqH_mass_FEM = [46.187, 122.90, 226.91, 357.35, 514.47];
freqL_hole_FEM = [49.752, 136.41, 257.08, 407.71, 585.07];
freqH_hole_FEM = [52.558, 149.50, 280.75, 454.81, 647.66];

%%%%%%%%%%
% figure %
%%%%%%%%%%
fig=figure(1);

gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];

%n
%freq0_mass
%freqL_mass
%freqH_mass
%freqL_hole
%freqH_hole

%%%%%%%%%%%%%
% subplot 1 %
%%%%%%%%%%%%%
subplot1 = subplot(2,1,1);

% p(1)=plot(n,freq0_mass,'-ko','DisplayName','perfect');
hold on;
p(3)=plot(n,freqL_mass,'--b','LineWidth',1.1,'DisplayName','1M, j = L (Analytic)');
p(4)=plot(n,freqL_mass_FEM,'bo','DisplayName','1M, j = L (FEM)');
p(5)=plot(n,freqH_mass,'--','LineWidth',1.1,'color',[0.58 0 0.83],'DisplayName','1M, j = H (Analytic)');
p(6)=plot(n,freqH_mass_FEM,'d','color',[0.58 0 0.83],'DisplayName','1M, j = H (FEM)');
p(7)=plot(n,freqL_hole,':r','LineWidth',1.5,'DisplayName','1H, j = L (Analytic)');
p(8)=plot(n,freqL_hole_FEM,'rs','DisplayName','1H, j = L (FEM)');
p(9)=plot(n,freqH_hole,':','LineWidth',1.5,'color',[1.0 0.5 0.0],'DisplayName','1H, j = H (Analytic)');
p(10)=plot(n,freqH_hole_FEM,'^','color',[1.0 0.5 0.0],'DisplayName','1H, j = H (FEM)');
p(1)=plot(n,freq0_mass,'-k','LineWidth',.7,'DisplayName','perfect (Analytic)');
p(2)=plot(n,freq0_FEM,'k+','DisplayName','perfect (FEM)');

set(gca,'XTick',2:6);
set(gca,'XTickLabel',{'','','','',''});

xlim([2.-0.001 6.+0.001]);

%xl1=xlabel('Erase');
yl1=ylabel('Frequency [Hz]','Interpreter','Latex');

legend(p(1:10),'Location','NorthWest'); legend boxoff;

%%%%%%%%%%%%%
% subplot 2 %
%%%%%%%%%%%%%
subplot2 = subplot(2,1,2);

hold on;
p(11)=plot(n,freqL_mass./freq0_mass,'--b','LineWidth',1.1,'DisplayName','1M, j = L (Analytic)');
p(12)=plot(n,freqL_mass_FEM./freq0_FEM,'bo','DisplayName','1M, j = L (FEM)');
p(13)=plot(n,freqH_mass./freq0_mass,'--','LineWidth',1.1,'color',[0.58 0 0.83],'DisplayName','1M, j = H (Analytic)');
p(14)=plot(n,freqH_mass_FEM./freq0_FEM,'d','color',[0.58 0 0.83],'DisplayName','1M, j = H (FEM)');
p(15)=plot(n,freqL_hole./freq0_mass,':r','LineWidth',1.5,'DisplayName','1H, j = L (Analytic)');
p(16)=plot(n,freqL_hole_FEM./freq0_FEM,'rs','DisplayName','1H, j = L (FEM)');
p(17)=plot(n,freqH_hole./freq0_mass,':','LineWidth',1.5,'color',[1.0 0.5 0.0],'DisplayName','1H, j = H (Analytic)');
p(18)=plot(n,freqH_hole_FEM./freq0_FEM,'^','color',[1.0 0.5 0.0],'DisplayName','1H, j = H (FEM)');

set(gca,'XTick',2:6);
set(gca,'XTickLabel',{'2','3','4','5','6'});

set(gca,'YTick',0.8:0.2:1.4);
set(gca,'YTickLabel',{'0.8','1.0','1.2',''});

xlim([2.-0.001 6.+0.001]);

xl2=xlabel('mode number','Interpreter','Latex');
yl2=ylabel('Frequency Ratio ($r$)','Interpreter','Latex');

%%%%%%%%%%%%%%%%%%%%%%
% subplot properties %
%%%%%%%%%%%%%%%%%%%%%%
set(subplot1, 'position', [0.12,0.37,.85,.6]);
set(subplot2,'position',[0.12,0.1,.85,.25]);

subplot1;
subplot2;

grid(subplot1,'on');
grid(subplot2,'on');

subplot1.FontSize = 10;
subplot2.FontSize = 10;
%xl1.FontSize = 1;
%yl1.FontSize = 10;
%xl2.FontSize = 10;
%yl2.FontSize = 10;

%%%%%%%%%%%%%
% save plot %
%%%%%%%%%%%%%
saveas(fig,['1imp_modeN.png']);
saveas(fig,['1imp_modeN.pdf']);


fprintf('\n\n');

%%% clear figure %%%
clf(figure(1));

fprintf('# ----- END Calculation ----- # \n');
quit;

