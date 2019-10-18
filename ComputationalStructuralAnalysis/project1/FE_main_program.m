%%  2D Plane Stress Element for HW#6

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 17 2017
%  Computational Structural Analysis Spring 2017

%--- V A R I A B L E   D E F I N I T I O N---------------------------------
%   'title': title of project
%   'npoin': number of nodes
%%%   'nelem': number of elements
%%%   'nnode': number of nodes per element
%   'nmate': number of materials
%   'nsect': number of sections
%%%   'ndofn': number of DOFs per node
%   'nevab': size of elemental stiffness matrix
%   'lnods': [ node1 node2 node3 node4;
%            node1 node2 node3 node4;
%            ... ];  each row is for each element
%   Note size of 'lnods' is (number of nodes x number of nodes per element)
%--- V A R I A B L E   D E F I N I T I O N---------------------------------

clear all; clc;

global nelem nnode ndofn

%-----------------------------------------
%  Read FE data from input data file
%-----------------------------------------

% open input file that contains FE data (nodal coord. & element
% connectivity & sectional properties & material properties)
ifid = fopen('FEdata.inp','rt');

% scan input file
ScanInput(ifid);

%-----------------------------------------
%  Compute Element Stiffness Matrix
%-----------------------------------------

% compute nevab: size of elemental stiffness matrix
nevab = nnode*ndofn; % nnode: number of nodes per element, ndofn: number of DOFs per node

% initialize element stiffness matrix
estif = zeros(nevab,nevab);

% determine gauss quadrature rule to use
ngaus = 2; % two point rule

% set Gauss Quadrature Constants
SetGaussQ(ngaus);

for ielem=1:nelem
    % clear element stiffness matrix
    estif = zeros(nevab,nevab);

    % calculate linear membrane stiffness matrix
    [estif] = GetElemStiff(ielem,ngaus,estif)
end

fclose('all');
quit;
