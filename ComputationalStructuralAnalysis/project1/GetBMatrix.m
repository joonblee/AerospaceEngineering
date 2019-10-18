function [Bmat] = GetBMatrix(Jmat,SFD,Bmat)
%% compute B martix

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 17 2017
%  Computational Structural Analysis Spring 2017

%	A = invJ*sfd;			[ N1,x N2,x N3,x N4,x ]	-> 4-node element
%							[ N1,y N2,y N3,y N4,y ]

%   B_m matrix		[N1,x 0     ..... N4,x 0    ]	-> 4node element
%					[0    N1,y  ..... 0    N4,y ]
%					[N1,y N1,x  ..... N4,y N4,x ]

% compute inverse of Jacobian matrix
% inv(Jmat)

% compute invJ * SFD
Np = inv(Jmat)*SFD; % = N'

% construct B matrix
Bmat=[Np(1,1) 0       Np(1,2) 0       Np(1,3) 0       Np(1,4) 0;
      0       Np(2,1) 0       Np(2,2) 0       Np(2,3) 0       Np(2,4);
	  Np(2,1) Np(1,1) Np(2,2) Np(1,2) Np(2,3) Np(1,3) Np(2,4) Np(1,4)];

return
