%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.3 Mode shape of the imperfect shell                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('--------------------- Start Code -----------------------\n');

C_i = zeros(length(phi_i),1);
denominator = 0;
numerator   = 0;

for i = 1:length(theta_i)
  fprintf(['i =' i '\n']);
  C_i(i) = m_i(i)*tan(phi_i(i)/2)^(2*n)*(n+cos(phi_i(i)))^2
  numerator   = numerator + C_i(i)*sin(2*n*theta_i(i))
  denominator = denominator + C_i(i)*cos(2*n*theta_i(i))
end

%%% shift angle %%%
j=[0 2];
zetaL = 1/(2*n) * (atan(numerator/denominator)+j*pi) %( theta_i + pi ); % j = 0 = L
zetaH = 1/(2*n) * (atan(numerator/denominator)+(j+1)*pi) %theta_i; % j = 1 = H

%%% displacement in alpha direction %%%
wL = zeros(2,length(theta));
wH = zeros(2,length(theta));

for i=1:2
  for j=1:length(theta)
    wL(i,j) = a+D*a*tan(phi/2)^n*(n+cos(phi))*sin(n*(theta(j)-zetaL(i)));
    wH(i,j) = a+D*a*tan(phi/2)^n*(n+cos(phi))*sin(n*(theta(j)-zetaH(i)));
  end
end

%%% plot %%%
ftitle = ['n = ' num2str(n)]; % ['Mode shape of the imp shell (mode ' num2str(n) ')'];
% fname  = '';
fsign  = '';
fphi   = '';
ftheta = '';

% char(10)
for i = 1:length(theta_i)
  % ftitle = [ftitle char(10) '(m,\phi,\theta' sprintf(')_%d = ',i) '(' num2str(m_i(i)) ',' num2str(rad2deg(phi_i(i))) '^{\circ},' num2str(rad2deg(theta_i(i))) '^{\circ})'];
  % fname = [fname 'm' num2str(m_i(i)) 'p' num2str(rad2deg(phi_i(i))) 't' num2str(rad2deg(theta_i(i)))];
  fsign = [fsign SIGN_i(i)];
  fphi  = [fphi '-' num2str(rad2deg(phi_i(i)))];
  ftheta = [ftheta '-' num2str(rad2deg(theta_i(i)))];
end

figure(1);
plot_radius=polarplot(theta,a*ones(1,length(theta)),'color',[0 0.4470 0.7410]);
plot_radius.LineWidth = 0.1;

% polar axes
pax = gca;
%pax.ThetaAxisUnits = 'radians';
pax.RMinorGrid = 'on';

pax.ThetaTickLabel = [{'\theta = 0^{\circ}','\theta = 30^{\circ}','\theta = 60^{\circ}','\theta = 90^{\circ}','\theta = 120^{\circ}','\theta = 150^{\circ}','\theta = 180^{\circ}','\theta = 210^{\circ}','\theta = 240^{\circ}','\theta = 270^{\circ}','\theta = 300^{\circ}','\theta = 330^{\circ}'}] 

for i=1:2
  figure(1);
  hold on;
  plot_wL(i)=polarplot(theta,wL(i,:),'b');
  plot_wH(i)=polarplot(theta,wH(i,:),'r--');
end

% locations of attached masses/holes
location_mass = [];
location_hole = [];
theta_mass = [];
theta_hole = [];

for i=1:length(m_i)
  figure(1);
  hold on;
  if m_i(i)>0
    theta_mass=[theta_mass,theta_i(i)];
    location_mass=[location_mass,a*sin(phi_i(i))];
  else
    theta_hole=[theta_hole,theta_i(i)];
    location_hole=[location_hole,a*sin(phi_i(i))];
  end
end

figure(1);
hold on;
plot_mass=polarplot(theta_mass,location_mass,'o','color',[0.6350 0.0780 0.1840]);
if length(location_hole)>0
  plot_hole=polarplot(theta_hole,location_hole,'x','color',[0.9290 0.6940 0.1250]);
end

figure(1);
if length(location_hole)==0
  leg=legend([plot_mass plot_wL(1) plot_wH(1)], 'attached mass', 'j = L', 'j = H' ,'Position',[0.1 0.1 0.1 0.1]);
elseif length(location_mass)==0
  leg=legend([plot_hole plot_wL(1) plot_wH(1)], 'hole', 'j = L', 'j = H' ,'Position',[0.1 0.1 0.1 0.1]);
else
  leg=legend([plot_mass plot_hole plot_wL(1) plot_wH(1)], 'attached mass', 'hole', 'j = L', 'j = H' ,'Position',[0.1 0.1 0.1 0.1]);
end
%legend boxoff;


title( ftitle );

saveas(figure(1),['figures/fig4_n',num2str(n),fsign,'_p',fphi,'_t',ftheta,'.png']);

fprintf('\n\n');

%%% clear figure %%%
clf(figure(1));

fprintf('# ----- END Calculation ----- # \n');
quit;

