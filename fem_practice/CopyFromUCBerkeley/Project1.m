%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ghaida El-Saied
%Project 1
%Finite Element Analysis
%University of California at Berkeley

% Force is given
% A program to solve the fh = -(k^2*sin((pi*k*xh)/L))/A - (2*xh)/A;
% Boundary value problem

% Conditions
% global stiffness = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Boundary conditions
%L=length, A = amplitude, u0 and uL are first and second boundary
%conditions respectively

A = 0.2; L= 1; u0 = 0; uL = 1; 

%BVP
%fh = -(k^2*sin((pi*k*xh)/L))/A - (2*xh)/A;
%true solution
%u_t = (x/L) + (x/(3.*A)).*(x.^2 - L.^2) - (L.^2/(A.*pi.^2)).*sin((pi.*k.*x)/L); 


%k and ne (number of elements) values:
k = [1 2 4 8 16 32]; % test value of 'wave num/2' (wave number = 1/wave length)
ne = [2 4 8 16 32 64 128 256 512 1024 2048]; % # of nodes
test_array = [1 1 1 1 1 1 1 1 1 1 1];


%setting up array for plotting
count = 1; 
label = cell(1,count);
bar_x = zeros(1,6);
bar_y = zeros(1,6); 
i=1;

%plot style
plot_style = ["-","-bo","-rs","-m^","-cx","--g+","--y*","k."];


%displaying figure for k values - is outputted below
for ii = 1:length(k)
  figure(ii);     
  
  for j = 1:length(ne) % Ne values to loop around
    %Global Matrices: stiffness = K, Force = F, Displacement = D
    K = zeros(ne(j) + 1, ne(j) + 1); % global stiffness matrix
    F = zeros(ne(j) + 1, 1); 
    D = zeros(ne(j) + 1, 1);

    % create a mesh
    he = 1/ne(j); % length of each elements
    xh = 0:he:L; % location of each nodes
    
    fh = -(k(ii)^2*sin((pi*k(ii)*xh)/L))/A - (2*xh)/A;

    e =1;
    while e <= ne(j) %loops over the ne matrix
      K_local = 1/he * [1 -1; -1 1]; % local stiffness matrix (stiffness = E*A/L, In this prob., total stiff of the beam is 1. local stiffness = 1/L)

      % discrete forces at each node:
      f_left = fh(e); f_right = fh(e+1);

      % force for each element: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      f_e = he/6 * [2 1;1 2] * [f_left; f_right]; 

      %creating force vectors
      if e == 1 
        f_e = f_e - u0 * K_local(:,1);
      end
      if e == ne(j)
        f_e = f_e - uL * K_local(:,2); 
      end

      %combining stiffness and force elemennts
      K(e:e+1, e: e+1) = K(e:e+1, e: e+1) + K_local;
      F(e:e+1) = F(e:e+1) + f_e;
      e = e + 1;
    end
        
    %Displacement
    D(1) = 0;
    D(ne(j) + 1) = 1;
    D(2:ne(j)) = K(2:ne(j), 2:ne(j))\F(2:ne(j));  
   
    plot(xh, D, plot_style(j),'Linewidth',1)
    hold on;

    %calculating error for each # of elements and k value  
    x = 0:L/1000:L;
    u_h = interp1(xh,D,x);
    u_t = (x/L) + (x/(3.*A)).*(x.^2 - L.^2) - (L.^2/(A.*pi.^2)).*sin((pi.*k(ii).*x)/L);
    diff = (u_t - u_h).^2; 
	error = trapz(x,diff);
    label{1,count} = ['Ne = ' num2str(ne(j)) ', '  'error = ' num2str(error)];
    if error <= 0.01
      % tolerance = error;
      bar_x(1,i) = k(ii);
      bar_y(1,i) = ne(j); 
      count = 1;
      i = i+1;
      break;
    end
    count = count +1;
  end

  %plot for true solution:
  
  hold on;
  plot(x,u_t,':','Linewidth',1)
  legend([label,'u^t true solution']);
  legend('Location', 'NorthWestOutside');
  xlabel('x'); ylabel('displacement');

  saveas(figure(ii),['figures_Project1/fig' num2str(ii) '.png']);
end

figure(length(k)+1);;
cat = categorical(bar_x);
bar(cat,bar_y); 
grid on;
xlabel('K');ylabel('number of elements ne');
title('K vs. ne');

saveas(figure(length(k)+1),['figures_Project1/fig' num2str(length(k)+1) '.png']);

fprintf('# ===== SUCCESS to run "Project1.m" ===== #\n\n');

quit;
