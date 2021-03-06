%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trimming of imperfect hemispherical shell including point mass distributions %
% Fig.3 Mode shape of the imperfect shell                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% variables %%%
phi   = pi/2;
theta = [0:359]/180*pi; 

%%% constants %%%
a=10; % radius
D=.1; % mode shape amplitude


%%% shift angle %%%
zetaL = 1/(2*n) * ( theta_i + pi ); % j = 1 = L
zetaH = 1/(2*n) * theta_i; % j = 2 = H

%%% displacement in alpha direction %%%
wL = a+D*a*tan(phi/2)^n*(n+cos(phi))*sin(n*(theta-zetaL))
wH = a+D*a*tan(phi/2)^n*(n+cos(phi))*sin(n*(theta-zetaH))

%%% plot %%%
figure(1);
title('Mode shape of the imperfect shell');

figure(1);
plot_shape(1) = polarplot(theta,wL);

% polar axes
pax = gca;
%pax.ThetaAxisUnits = 'radians';
pax.RMinorGrid = 'on';

figure(1);
hold on;
plot_shape(2) = polarplot(theta,wH,'r--');

figure(1);
legend(['mode # =',num2str(n),', \phi =',num2str(rad2deg(phi_i)),'^{\circ} , \theta =',num2str(rad2deg(theta_i)),'^{\circ} , j = L'],['mode # =',num2str(n),', \phi =',num2str(rad2deg(phi_i)),'^{\circ} , \theta =',num2str(rad2deg(theta_i)),'^{\circ} , j = H']);
legend boxoff;

saveas(figure(1),['figures/fig3_n',num2str(n),'phi',num2str(rad2deg(phi_i)),'theta',num2str(rad2deg(theta_i)),'.png']);

fprintf('\n\n');

%%% clear figure %%%
clf(figure(1));

fprintf('# ----- END Calculation ----- # \n');
quit;

