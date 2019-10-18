%%  2D Plane Stress Element for HW#6

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

%--- V A R I A B L E   D E F I N I T I O N---------------------------------
%   'title': title of project
%   'npoin': number of nodes
%   'nelem': number of elements
%   'nnode': number of nodes per element
%   'nmate': number of materials
%   'nsect': number of sections
%   'ndofn': number of DOFs per node
%   'nevab': size of elemental stiffness matrix
%   'nvfix': number of restrained nodes
%   'lnods': [ node1 node2 node3 node4;
%            node1 node2 node3 node4;
%            ... ];  each row is for each element
%   Note size of 'lnods' is (number of nodes x number of nodes per element)
%   'IDArray': [ ndofn x npoin ] matrix that stores equation number
%   'iffix': vector of size [1 x npoin*ndofn],which store 0 for active DOF 
%            and stores 1 for inactive DOF
%            // iffix = zeros(1,npoin*ndofn);
%   'neq': Total number of equations, that is, total number of active DOFs
%--- V A R I A B L E   D E F I N I T I O N---------------------------------

clear all;
% clc;

global nelem nnode ndofn nevab neq
global lnods
global IDArray

%-----------------------------------------
%  Read FE data from input data file
%-----------------------------------------

% open input file that contains FE data (nodal coord. & element
% connectivity & sectional properties & material properties)
ifid = fopen('FEdata.inp','rt');

% scan input file
ScanInput(ifid);

% Get IDArray
IDArray = GetIDArray();

%-----------------------------------------
%  Compute Stiffness Matrix
%-----------------------------------------

% compute size of elemental stiffness matrix
nevab = nnode*ndofn;

% initialize element stiffness matrix
estif = zeros(nevab,nevab);
xGK = zeros(neq,neq); % 

% determine gauss quadrature rule to use
ngaus = 2; % two point rule

% set Gauss Quadrature Constants
SetGaussQ(ngaus);

for ielem=1:nelem
    
    % clear element stiffness matrix
    estif = zeros(nevab,nevab);

    % calculate linear membrane stiffness matrix
    [estif] = GetElemStiff(ielem,ngaus,estif);
    
    % when using IDArray
    for inode=1:nnode
        nod = lnods(ielem,inode); % node number
    	for jdofn=1:ndofn
    		pos =(inode-1)*ndofn + jdofn; % numbering all dofn in this element
    		KK(pos) = IDArray(jdofn, nod); 
			% 'IDArray': [ ndofn x npoin ] matrix that stores equation number
			% IDArray gives (1) 0 for Non-active DOF
			%               (2) numbering for all Active DOFs in the whole structure
			%                   from 0 to nevab
        end
    end

    % Assemble element stiffness matrix
    for ievab=1:nevab
        if( KK(ievab) <= 0 )
            continue;
        end
        I=KK(ievab); % ID numbering
        for jevab=1:nevab
            J = KK(jevab);
            if( J < I )
                continue;
            end
            xGK(I,J) = xGK(I,J) + estif(ievab,jevab); % Global Stiffness Matrix
        end
    end
    
end

% fprintf('\n +++ global stiff mat +++ \n');
% xGK

% make it symmetric matrix
for ieq=1:neq
    for jeq=ieq:neq
        xGK(jeq,ieq) = xGK(ieq,jeq);
    end
end

%-----------------------------------------------
%  Compute Load Vector Matrix
%-----------------------------------------------
[GL] = GetDLoad(ngaus,IDArray);

%-----------------------------------------------
%  Solve for Nodal Displacements
%-----------------------------------------------
Disp = xGK\GL

%-----------------------------------------------
%  Compute stresses and strains at Gauss Points
%-----------------------------------------------
[Disp,Strs,Strn] = GetStrsStrn(ielem,ngaus,xGK,GL);

%fclose('all');
quit();
