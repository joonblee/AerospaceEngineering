function [coord2D] = GetCoord(ielem,coord2D)
%% get elemental coordinate matrix

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

global lnods nnode
global COORD

for inode=1:nnode
    node = lnods(ielem,inode);
    coord2D(inode,1) = COORD(node,2);
    coord2D(inode,2) = COORD(node,3);
end


return




