%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.4 Frequencies with respect to the ratio of imperfection                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%% constants %%%
n=2;

E=210E9; rho=7800; mu=0.3;

r= [0:.01:.1];

freq0_mass = [];
freq0_hole = [];
freqL_mass = [];
freqH_mass = [];
freqL_hole = [];
freqH_hole = [];

for pm = [1 -1]
  M = 5; %5./(1+pm*r)

  phi0 = pi/2;

  hOa = 0.01; % a/h
  a=((M*3/pi/rho)/(2+cos(phi0))/(1-cos(phi0))^2./(3*hOa+(hOa).^3/4)).^(1/3);
  h=a*hOa; % width

  phi_i   = pi/2;
  theta_i =    0;

  %%% Run iteration %%%
  % attached mass
  m_i = pm*r.*M;

  % variables
  syms phi;

  % shift angle
  zetaL = 1/(2*n) * theta_i; % j = 0 = L
  zetaH = 1/(2*n) * (theta_i + pi ); % j = 1 = H

  % integrations
  U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,phi0);
  K_phi = double(int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,phi0));

  K_iL = tan(phi_i/2)^(2*n) * (sin(phi_i)^2 + (n+cos(phi_i))^2 * sin(n*(theta_i-zetaL))^2 );
  K_iH = tan(phi_i/2)^(2*n) * (sin(phi_i)^2 + (n+cos(phi_i))^2 * sin(n*(theta_i-zetaH))^2 );


  % epsilon
  epsilonK_L = K_iL/K_phi*m_i./(pi*rho*a.^2.*h);
  epsilonK_H = K_iH/K_phi*m_i./(pi*rho*a.^2.*h);

  % eigen frequency (angular freq, natural frequency)
  omega0 = double(( n^2*(n^2-1)^2 * E*h.^2 ./ (3*(1+mu)*rho*a.^4) * U_phi/K_phi ).^.5);
  omegaL = double(omega0./(1+epsilonK_L).^.5);
  omegaH = double(omega0./(1+epsilonK_H).^.5);

  % linear frequency
  freq0 = double(omega0 / (2*pi));
  freqL = double(omegaL / (2*pi));
  freqH = double(omegaH / (2*pi));

  if pm == 1
    freq0_mass = freq0(1) * ones(1,length(r));
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
Xaxis_FEM = [0.02,0.04,0.06,0.08,0.1];
freq0_mass_FEM = 49.164;
freqL_mass_FEM = [48.568,48.052,47.601,47.203,46.850];
freqH_mass_FEM = [46.270,43.921,41.986,40.371,39.005];
freq0_hole_FEM = 49.147;
freqL_hole_FEM = [49.860,50.682,51.668,52.870,54.362];
freqH_hole_FEM = [52.663,56.848,61.525,66.286,70.721];


%%%%%%%%%%
% Figure %
%%%%%%%%%%

fig=figure(1);

gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];

grid on;

hold on;
p(1)=plot(r,freqL_mass./freq0_mass,'b','LineStyle','--','LineWidth',.5,'DisplayName','1M, j = L (Analytic)');
p(2)=plot(Xaxis_FEM,freqL_mass_FEM/freq0_mass_FEM,'bo','MarkerSize',3.,'DisplayName','1M, j = L (FEM)');
p(3)=plot(r,freqH_mass./freq0_mass,'color',[0.58 0 0.83],'LineStyle','--','LineWidth',.5,'DisplayName','1M, j = H (Analytic)');
p(4)=plot(Xaxis_FEM,freqH_mass_FEM/freq0_mass_FEM,'d','MarkerSize',3.,'color',[0.58 0 0.83],'DisplayName','1M, j = H (FEM)');
p(5)=plot(r,freqL_hole./freq0_mass,'r','LineStyle',':','LineWidth',1.1,'DisplayName','1H, j = L (Analytic)');
p(6)=plot(Xaxis_FEM,freqL_hole_FEM/freq0_hole_FEM,'rs','MarkerSize',3.,'DisplayName','1H, j = L (FEM)');
p(7)=plot(r,freqH_hole./freq0_mass,'color',[1.0 0.5 0.0],'LineStyle',':','LineWidth',1.1,'DisplayName','1H, j = H (Analytic)');
p(8)=plot(Xaxis_FEM,freqH_hole_FEM/freq0_hole_FEM,'^','MarkerSize',3.,'color',[1.0 0.5 0.0],'DisplayName','1H, j = H (FEM)');


xlabel('$m_{imp}$ / $m_{shell}$','Interpreter','Latex');
ylabel('Frequency Ratio ($r$)','Interpreter','Latex');

legend(p(1:8),'Location','northwest'); legend boxoff;

saveas(fig,'1imp_massRatio.png');
saveas(fig,'1imp_massRatio.pdf');

%%% End of the code %%%
fprintf('# ----- END Calculation ----- # \n');

quit;
