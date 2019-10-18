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

hOa=[0.01:0.01:0.1] % a/h
a=((15/(2*pi*rho))./(3*hOa+(hOa).^3/4)).^(1/3) % radius
% 1/2 rho [ 4/3 pi (a+h/2)^3 - 4/3 pi (a-h/2)^3 ] = 5 kg
h=a.*hOa % width

mass=4/3*pi*rho.*((a+h/2).^3-(a-h/2).^3)/2

%%% integrations %%%
U_phi = int( tan(phi/2)^(2*n) / sin(phi)^3 ,phi,0,pi/2);
K_phi = int( tan(phi/2)^(2*n) * ((n+cos(phi))^2+2*sin(phi)^2) * sin(phi) ,phi,0,pi/2);

%%% natural frequency of the perfect HS (omega0) %%%
omega0 = double(( n^2*(n^2-1)^2 * E*h.^2 ./ (3*(1+mu)*rho*a.^4) * U_phi/K_phi ) .^ .5) % angular freq
NatFreq0 = double(omega0 / (2*pi)) % linear freq

%%% figure %%%
figure(1);
grid on;
% title('Eigenfrequencies of the perfect shell'); 
xlabel('h/a'); 
ylabel('Frequency [Hz]');  %('Natural Frequency (=\omega_0/2\pi) [Hz]');

%%% plot (theory) %%%
figure(1);
hold on;
plot1(1)=plot(hOa,NatFreq0,'-','DisplayName','Analytical');

%%% plot (FEM) %%%
COMSOLx=hOa;
COMSOLy=[49.849;120.42;205.77;300.89;395.49;504.42;606.35;722.68;836.59;954.85];

figure(1);
hold on;
plot1(2)=plot(COMSOLx,COMSOLy,'x','DisplayName','FEM (COMSOL)');

xlim([0.01 0.1])

%%% save figure %%%
figure(1); legend(plot1(1:2),'Location','best'); legend boxoff;

saveas(figure(1),'figures/fig2.png');


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
