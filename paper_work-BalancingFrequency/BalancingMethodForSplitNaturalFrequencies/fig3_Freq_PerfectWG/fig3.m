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

M = 5 % mass of the structure
hOa=[0.01:0.01:0.1] % a/h
a=((M*3/pi/rho)/(2+cos(phi0))/(1-cos(phi0))^2./(3*hOa+(hOa).^3/4)).^(1/3) % radius
% 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg

h=a.*hOa % width

mass=pi*rho/3*((a+h/2).^3-(a-h/2).^3)*(2+cos(phi0))*(1-cos(phi0))^2
% https://en.wikipedia.org/wiki/Spherical_cap
% It should be same with 'M'.

%%% integrations %%%
U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,phi0);
K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,phi0);

%%% natural frequency of the perfect HS (omega0) %%%
omega0 = double(( n^2*(n^2-1)^2 * E*h.^2 ./ (3*(1+mu)*rho*a.^4) * U_phi/K_phi ) .^ .5) % angular freq
NatFreq0 = double(omega0 / (2*pi)) % linear freq

  if phi0 == pi/4
    freq_45 = NatFreq0;
  elseif phi0 == pi/2
	freq_90 = NatFreq0;
  else
	freq_120 = NatFreq0;
  end

end


%%% figure %%%
fig=figure(1);

gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];

grid on;
% title('Eigenfrequencies of the perfect shell'); 
xlabel('h/a'); 
ylabel('Frequency [Hz]');  %('Natural Frequency (=\omega_0/2\pi) [Hz]');

%%% plot (FEM) %%%
COMSOLx=[0.01;0.04;0.07;0.1];%hOa;
COMSOLy_45=[57.534;354.77;727.55;1143.0];
COMSOLy_90=[49.849;300.89;606.35;954.85];
COMSOLy_120=[69.095;414.52;829.55;1270.4];

%%% plot (theory) %%%
%figure(1);
fig;
hold on;
plot1(1)=plot(hOa,freq_45,'-','color',[0 0 1],'LineWidth',1.1,'DisplayName','\phi_0=\pi/4 (Analytic)');

%figure(1);
fig;
hold on;
plot1(2)=plot(COMSOLx,COMSOLy_45,'x','color',[0 0.4470 0.7410],'DisplayName','\phi_0=\pi/4 (FEM)');

%figure(1);
fig;
hold on;
plot1(3)=plot(hOa,freq_90,'--','color',[1 0 0],'LineWidth',1.1,'DisplayName','\phi_0=\pi/2 (Analytic)');

%figure(1);
fig;
hold on;
plot1(4)=plot(COMSOLx,COMSOLy_90,'o','color',[0.8500 0.3250 0.0980],'DisplayName','\phi_0=\pi/2 (FEM)');

%figure(1);
fig;
hold on;
plot1(5)=plot(hOa,freq_120,':','color',[1 0 1],'LineWidth',1.1,'DisplayName','\phi_0=2\pi/3 (Analytic)');

%figure(1);
fig;
hold on;
plot1(6)=plot(COMSOLx,COMSOLy_120,'^','color',[0.4940 0.1840 0.5560],'DisplayName','\phi_0=2\pi/3 (FEM)');

xlim([0.01 0.1])

%%% save figure %%%
%figure(1); 
fig; legend(plot1(1:6),'Location','northwest'); legend boxoff;

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
