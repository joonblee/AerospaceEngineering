function [Bmat] = GetBMatrix(Jmat,SFD,Bmat)
%% compute B martix

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

%	A = invJ*sfd;			[ N1,x N2,x N3,x N4,x ]	-> 4-node element
%							[ N1,y N2,y N3,y N4,y ]

%   B_m matrix		[N1,x 0     ..... N4,x 0    ]	-> 4node element
%					[0    N1,y  ..... 0    N4,y ]
%					[N1,y N1,x  ..... N4,y N4,x ]

% compute inverse of Jacobian matrix
invJ = inv(Jmat);

% compute invJ * SFD
A = invJ*SFD;

% construct B matrix
Bmat(1,1) = A(1,1);
Bmat(2,2) = A(2,1);
Bmat(3,1) = A(2,1);
Bmat(3,2) = A(1,1);

Bmat(1,3) = A(1,2);
Bmat(2,4) = A(2,2);
Bmat(3,3) = A(2,2);
Bmat(3,4) = A(1,2);

Bmat(1,5) = A(1,3);
Bmat(2,6) = A(2,3);
Bmat(3,5) = A(2,3);
Bmat(3,6) = A(1,3);

Bmat(1,7) = A(1,4);
Bmat(2,8) = A(2,4);
Bmat(3,7) = A(2,4);
Bmat(3,8) = A(1,4);

return