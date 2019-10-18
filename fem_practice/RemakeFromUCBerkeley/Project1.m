% Force is given
% A program to solve the fi = -(k^2*sin((pi*k*xi)/L))/A - (2*xi)/A;
% Boundary value problem

% Conditions
% global stiffness = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Boundary conditions
%L=length, A = amplitude, u0 and uL are first and second boundary
%conditions respectively

A = 0.2; L= 1; u0 = 0; uL = 1; 

%BVP
%fi = -(k^2*sin((pi*k*xi)/L))/A - (2*xi)/A;
%true solution
%u_t = (x/L) + (x/(3.*A)).*(x.^2 - L.^2) - (L.^2/(A.*pi.^2)).*sin((pi.*k.*x)/L); 


%k and ne (number of elements) values:
k = 8; % test value of 'wave num/2' (wave number = 1/wave length)
Ne = [2 4 8 16 32 64 128 256 512 1024 2048]; % # of nodes
test_array = [1 1 1 1 1 1 1 1 1 1 1];


%setting up array for plotting
count = 1; 
label = cell(1,count);
i=1;

%plot style
plot_style = ["-","-bo","-rs","-m^","-cx","--g+","--y*","k."];


%displaying figure for k values - is outputted below
figure(1);     
  
for j = 1:length(Ne) % Ne values to loop around
  K_global = zeros(Ne(j) + 1, Ne(j) + 1); % global stiffness matrix
  F_global = zeros(Ne(j) + 1, 1);         % global force vector
  D = zeros(Ne(j) + 1, 1);                % displacement
  D_ = zeros(Ne(j)+1,1);

  % create a mesh
  L_e = 1/Ne(j);                             % length of each element
  k_e = 1/L_e;                               % stiffness of each element
  xi = 0:L_e:L;                              % location at each node
  fi = -(k^2*sin((pi*k*xi)/L))/A - (2*xi)/A; % descrete forces at each node

  for e = 1:Ne(j) %loops over the Ne (# of elements)
    K_local = k_e * [1 -1; -1 1]; % local stiffness matrix (stiffness = E*A/L, In this prob., total stiff of the beam is 1. local stiffness = 1/L)
    K_global(e:e+1, e:e+1) = K_global(e:e+1, e:e+1) + K_local;


    % force for each element: 
    f_e = L_e/6 * [2 1;1 2] * [fi(e); fi(e+1)]; 

    %creating force vectors
    if e == 1 
      f_e = f_e - u0 * K_local(:,1);
    elseif e == Ne(j)
      f_e = f_e - uL * K_local(:,2); 
    end

    %combining force elemennts
    F_global(e:e+1) = F_global(e:e+1) + f_e;
  end

  %Displacement
  D(1) = 0;
  D(Ne(j) + 1) = 1;

  fprintf('\n------------------------------------\n');
  j
  D(2:Ne(j)) = K_global(2:Ne(j), 2:Ne(j))\F_global(2:Ne(j))
  D_ = K_global\F_global


  plot(xi, D, plot_style(j),'Linewidth',1)
  hold on;

  %calculating error for each # of elements and k value  
  x = 0:L/1000:L;
  u_h = interp1(xi,D,x);
  u_t = (x/L) + (x/(3.*A)).*(x.^2 - L.^2) - (L.^2/(A.*pi.^2)).*sin((pi.*k.*x)/L);
  diff = (u_t - u_h).^2; 
  error = trapz(x,diff);
  label{1,count} = ['Ne = ' num2str(Ne(j)) ', '  'error = ' num2str(error)];
  if error <= 0.01
    % tolerance = error;
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
legend('Location', 'best');
legend boxoff;
xlabel('x'); ylabel('displacement');

saveas(figure(1),['figures_Project1/fig.png']);

fprintf('# ===== SUCCESS to run "Project1.m" ===== #\n\n');

quit;
