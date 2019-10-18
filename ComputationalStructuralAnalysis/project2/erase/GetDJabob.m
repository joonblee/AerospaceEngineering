function [detJ,Jmat] = GetDJabob(SFD,coord2D,detJ,Jmat)
%% Get Jacobian matrix and its determinant value

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

Jmat = SFD*coord2D;
detJ = det(Jmat);

return