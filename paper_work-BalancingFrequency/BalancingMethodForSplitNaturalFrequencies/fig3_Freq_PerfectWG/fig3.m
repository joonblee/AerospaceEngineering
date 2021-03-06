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

freq_45=[]; freq_90=[]; freq_120=[];

for phi0 = [pi/4 pi/2 2*pi/3]

  M = 5; % mass of the structure
  hOa=[0.01:0.01:0.1]; % a/h
  a=((M*3/pi/rho)/(2+cos(phi0))/(1-cos(phi0))^2./(3*hOa+(hOa).^3/4)).^(1/3); % radius
  % 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg

  h=a.*hOa; % width

  mass=pi*rho/3*((a+h/2).^3-(a-h/2).^3)*(2+cos(phi0))*(1-cos(phi0))^2;
  % https://en.wikipedia.org/wiki/Spherical_cap
  % It should be same with 'M'.

  %%% integrations %%%
  U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,phi0);
  K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,phi0);

  %%% natural frequency of the perfect HS (omega0) %%%
  omega0 = double(( n^2*(n^2-1)^2 * E*h.^2 ./ (3*(1+mu)*rho*a.^4) * U_phi/K_phi ) .^ .5); % angular freq
  NatFreq0 = double(omega0 / (2*pi)); % linear freq

  if phi0 == pi/4
    freq_45 = NatFreq0;
  elseif phi0 == pi/2
	freq_90 = NatFreq0;
  else
	freq_120 = NatFreq0;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%
% FEM (COMSOL) results %
%%%%%%%%%%%%%%%%%%%%%%%%
COMSOLx=[0.01;0.04;0.07;0.1];%hOa;
COMSOLy_45=[57.534;354.77;727.55;1143.0];
COMSOLy_90=[49.849;300.89;606.35;954.85];
COMSOLy_120=[69.095;414.52;829.55;1270.4];

ref_freq_45=[];
ref_freq_90=[];
ref_freq_120=[];

for it_COMSOLx = 1:size(COMSOLx)
  ref_COMSOLx = COMSOLx(it_COMSOLx);
  for it_hOa = 1:10
    if hOa(it_hOa) == ref_COMSOLx
      ref_freq_45=[ref_freq_45; freq_45(it_hOa)];
      ref_freq_90=[ref_freq_90; freq_90(it_hOa)];
      ref_freq_120=[ref_freq_120; freq_120(it_hOa)];
      break;
    end
  end
end

%%%%%%%%%%
% figure %
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
p(1)=plot(hOa,freq_45,'-r','LineWidth',.5,'DisplayName','\phi_0=\pi/4 (Analytic)');
p(2)=plot(COMSOLx,COMSOLy_45,'xr','MarkerSize',4.,'DisplayName','\phi_0=\pi/4 (FEM)');
p(3)=plot(hOa,freq_90,'--b','LineWidth',.5,'DisplayName','\phi_0=\pi/2 (Analytic)');
p(4)=plot(COMSOLx,COMSOLy_90,'ob','MarkerSize',3.,'DisplayName','\phi_0=\pi/2 (FEM)');
p(5)=plot(hOa,freq_120,':','color',[0.58 0 0.83],'LineWidth',1.1,'DisplayName','\phi_0=2\pi/3 (Analytic)');
p(6)=plot(COMSOLx,COMSOLy_120,'^','color',[0.58 0 0.83],'MarkerSize',3.,'DisplayName','\phi_0=2\pi/3 (FEM)');

set(gca,'XTick',0.01:0.01:0.1);
set(gca,'XTickLabel',{'','','','','','','','','',''});

set(gca, 'YScale', 'log')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')
set(gca,'YTick',[50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000]);
set(gca,'YTickLabel',{'50', '60', '', '80', '', '100', '200', '300', '400', '500', '600', '', '800', '', '1000'});

xlim([0.01 0.1]);

ylabel('Frequency [Hz]','Interpreter','Latex');
legend(p(1:6),'Location','southeast'); legend boxoff;

%%%%%%%%%%%%%
% subplot 2 %
%%%%%%%%%%%%%
subplot2 = subplot(2,1,2);

hold on;
p(7)=plot(COMSOLx,(ref_freq_45-COMSOLy_45)./ref_freq_45,'-rx','LineWidth',.5,'MarkerSize',4.);
p(8)=plot(COMSOLx,(ref_freq_90-COMSOLy_90)./ref_freq_90,'--ob','LineWidth',.5,'MarkerSize',3.);
p(9)=plot(COMSOLx,(ref_freq_120-COMSOLy_120)./ref_freq_120,':','LineWidth',1.1,'color',[0.58 0 0.83]);

dummy=plot(COMSOLx,(ref_freq_120-COMSOLy_120)./ref_freq_120,'^','MarkerSize',3.,'color',[0.58 0 0.83]);

set(gca,'XTick',0.01:0.01:0.1);
set(gca,'XTickLabel',{'0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.1'});

set(gca,'YTick',0:0.05:0.2);

xlim([0.01 0.1]);
ylim([0. 0.21]);

xlabel('h/a','Interpreter','Latex'); 
ylabel('$|\frac{\textrm{Analytic} - \textrm{FEM}}{\textrm{Analytic}}|$','Interpreter','Latex');
% ylabel('$|\frac{Analytic - FEM}{Analytic}|$','Interpreter','Latex');

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

%%%%%%%%%%%%%%%
% save figure %
%%%%%%%%%%%%%%%
saveas(fig,'fig3.png');
saveas(fig,'fig3.pdf');


quit();


%%% Print results and errors %%%
fprintf('\n### Print results and errors (=1-Analytical/FEM) ###\n\n');
for i=1:length(COMSOLx)
  Analytical=0; err=0;
  for j=1:length(NatFreq0)
    if COMSOLx(i)==hOa(j)
      Analytical=NatFreq0(j);
      break;
    end
  end
  err=(1-Analytical/COMSOLy(i))*100;
  fprintf('x = %f , Analytical = %f , FEM = %f , err = %f \n\n',COMSOLx(i),COMSOLy(i),Analytical,err);
end

%%% End of the code %%%
fprintf('# ----- END Calculation ----- # \n');
quit;
