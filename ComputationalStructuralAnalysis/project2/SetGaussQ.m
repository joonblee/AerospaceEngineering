function SetGaussQ(ngaus)
%% set Gauss quadrature constants

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

global posgp weigp

% initialize variables
posgp = zeros(ngaus*ngaus,2);
weigp = zeros(ngaus*ngaus,2);

if ngaus==2
    posgp(1,1) = -0.5773502691896260;	% sqrt(1/3)
    posgp(1,2) = -0.5773502691896260;
    posgp(2,1) = 0.5773502691896260;
    posgp(2,2) = -0.5773502691896260;
    posgp(3,1) = 0.5773502691896260;
    posgp(3,2) = 0.5773502691896260;
    posgp(4,1) = -0.5773502691896260;
    posgp(4,2) = 0.5773502691896260;

    weigp(1,1) = 1.0;
    weigp(1,2) = 1.0;
    weigp(2,1) = 1.0;
    weigp(2,2) = 1.0;
    weigp(3,1) = 1.0;
    weigp(3,2) = 1.0;
    weigp(4,1) = 1.0;
    weigp(4,2) = 1.0;
end

% you can try other gauss quadrature rules


return