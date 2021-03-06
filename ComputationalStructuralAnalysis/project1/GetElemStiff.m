function [estif] = GetElemStiff(ielem,ngaus,estif)
%% compute elemental stiffness matrix

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 17 2017
%  Computational Structural Analysis Spring 2017

global nevab nnode
global ELEM_DATA
global SECTION_DATA
global MATERIAL_DATA

% set variables
xipos=0.0;
etapos=0.0;
xiweight=0.0;
etaweight=0.0;
detJ =0.0;
SFD = zeros(2,nnode);
coord2D = zeros(nnode,2);
Bmat = zeros(3,nevab);
Dmat = zeros(3,3);
Jmat = zeros(2,2);

% get material and section properties for the element
no_sect = ELEM_DATA(ielem,6);
no_mate = ELEM_DATA(ielem,7);

E = MATERIAL_DATA(no_mate,2);
Nu = MATERIAL_DATA(no_mate,3);
thick = SECTION_DATA(no_sect,2);

% get material constitutive matrix
% Dmat = E/(1+Nu)/(1-2*Nu)*[1-Nu Nu   0;
%                          Nu   1-Nu 0;
%						  0    0    (1-2*Nu)/2];

Dmat = E*thick/(1-Nu^2)*[1  Nu 0;
                         Nu 1  0;
						 0  0  (1-Nu)/2];

% loop over integration points in the xi direction
for igaus=1:ngaus
    
    % loop over integration points in the eta direction
    for jgaus=1:ngaus
        
        % initialize arrays
        Bmat = zeros(3,nevab);
        SFD = zeros(2,nnode);
        coord2D = zeros(nnode,2);
        
        % get sampling point location and associated weight factor (GetGPV.m)
		[xipos,etapos,xiweight,etaweight] = GetGPV(ngaus,igaus,jgaus,xipos,etapos,xiweight,etaweight);

        % compute shape function derivatives (GetSFD.m)
        [SFD] = GetSFD(SFD,xipos,etapos);

        % get coordinate values of the element (GetCoord.m)
		[coord2D] = GetCoord(ielem,coord2D);

        % get Jacobian matrix and determinant of Jmat (GetDJacob.m)
		[detJ,Jmat] = GetDJabob(SFD,coord2D,detJ,Jmat);

        % get B natrix (GetBMatrix.m)
		[Bmat] = GetBMatrix(Jmat,SFD,Bmat);

        estif = estif + Bmat'*Dmat*Bmat*detJ*xiweight*etaweight;
    end
end


