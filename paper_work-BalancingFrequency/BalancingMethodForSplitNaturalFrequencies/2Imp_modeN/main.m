%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab independent run %
%%%%%%%%%%%%%%%%%%%%%%%%%%

m_i     = [0.05 0.05];
phi_i   = [pi/2 pi/2];
theta_i = [   0 pi/2];


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

freqPF = [];
freqL_1m = [];
freqH_1m = [];
freqL_1h = [];
freqH_1h = [];
freqL_2m = [];
freqH_2m = [];
freqL_1m1h = [];
freqH_1m1h = [];
freqL_2h = [];
freqH_2h = [];

%%%%%%%%%%%%%
% 2I system %
%%%%%%%%%%%%%
for SignN = [1 2 3] 

  if SignN == 2
    m_i(2) = -m_i(2);
  elseif SignN == 3
    m_i(1) = -m_i(1);
  end

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


  if SignN == 1
    freqPF = freq0;
    freqL_2m = freqL;
    freqH_2m = freqH;
  elseif SignN == 2
    freqL_1m1h = freqL;
    freqH_1m1h = freqH;
  elseif SignN == 3
    freqL_2h = freqL;
    freqH_2h = freqH;
  end

end


%%%%%%%%%%%%%
% 1I system %
%%%%%%%%%%%%%
m_i     = [0.1];
phi_i   = [pi/2];
theta_i = [0];

for pm = [1 -1]
  m_i = pm*0.1;

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
    freqL_1m = freqL;
    freqH_1m = freqH;
  else
    freqL_1h = freqL;
    freqH_1h = freqH;
  end

end


%%%%%%%%%%
% Figure %
%%%%%%%%%%
fig=figure(1);

gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];

%%%%%%%%%%%%%
% subplot 1 %
%%%%%%%%%%%%%
subplot1 = subplot(2,1,1);

hold on;
p(2)=plot(n,freqL_1m,'-.b','LineWidth',.5,'DisplayName','1M');
p(3)=plot(n,freqH_1m,'-.b','LineWidth',.5,'DisplayName','1M');
p(4)=plot(n,freqL_1h,'-.r','LineWidth',.5,'DisplayName','1H');
p(5)=plot(n,freqH_1h,'-.r','LineWidth',.5,'DisplayName','1H');

p(6)=plot(n,freqL_2m,'--','LineWidth',.5,'color',[0.58 0 0.83],'DisplayName','2M');
p(7)=plot(n,freqH_2m,'--','LineWidth',.5,'color',[0.58 0 0.83],'DisplayName','2M');
p(8)=plot(n,freqL_1m1h,':','LineWidth',1.1,'color',[1.0 0.5 0.0],'DisplayName','1M1H');
p(9)=plot(n,freqH_1m1h,':','LineWidth',1.1,'color',[1.0 0.5 0.0],'DisplayName','1M1H');
p(10)=plot(n,freqL_2h,'--','LineWidth',.5,'color',[0 0.7 0],'DisplayName','2H');
p(11)=plot(n,freqH_2h,'--','LineWidth',.5,'color',[0 0.7 0],'DisplayName','2H');

% ods^><x*ph

p(1)=plot(n,freqPF,'-+k','DisplayName','perfect');


% ylim([0 700]);

set(gca,'XTick',2:6);
set(gca,'XTickLabel',{'','','','',''});

set(gca, 'YScale', 'log')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')
set(gca,'YTick',[50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000]);
set(gca,'YTickLabel',{'50', '60', '', '80', '', '100', '200', '300', '400', '500', '600', '', '800', '', '1000'});

ylabel('Frequency [Hz]','Interpreter','Latex');

legend(p([1 2 4 6 8 10]),'Location','NorthWest'); legend boxoff;

%%%%%%%%%%%%%
% subplot 2 %
%%%%%%%%%%%%%
subplot2 = subplot(2,1,2);

hold on;

p(2)=plot(n,freqL_1m./freqPF,'-.b','LineWidth',.5);
p(3)=plot(n,freqH_1m./freqPF,'-.b','LineWidth',.5);
p(4)=plot(n,freqL_1h./freqPF,'-.r','LineWidth',.5);
p(5)=plot(n,freqH_1h./freqPF,'-.r','LineWidth',.5);

p(6)=plot(n,freqL_2m./freqPF,'--','LineWidth',.5,'color',[0.58 0 0.83]);
p(7)=plot(n,freqH_2m./freqPF,'--','LineWidth',.5,'color',[0.58 0 0.83]);
p(8)=plot(n,freqL_1m1h./freqPF,':','LineWidth',1.1,'color',[1.0 0.5 0.0],'DisplayName','1M1H');
p(9)=plot(n,freqH_1m1h./freqPF,':','LineWidth',1.1,'color',[1.0 0.5 0.0],'DisplayName','1M1H');
p(10)=plot(n,freqL_2h./freqPF,'--','LineWidth',.5,'color',[0 0.7 0]);
p(11)=plot(n,freqH_2h./freqPF,'--','LineWidth',.5,'color',[0 0.7 0]);

set(gca,'XTick',2:6);
set(gca,'XTickLabel',{'2','3','4','5','6'});

xlabel('mode number','Interpreter','Latex');
ylabel('Frequency Ratio ($r$)','Interpreter','Latex');

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

%%%%%%%%%%%%%
% save plot %
%%%%%%%%%%%%%
saveas(fig,['2imp_modeN.png']);
saveas(fig,['2imp_modeN.pdf']);

fprintf('\n\n');

%%% clear figure %%%
clf(fig);

fprintf('# ----- END Calculation ----- # \n');
quit;

