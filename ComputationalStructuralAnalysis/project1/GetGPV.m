function [xipos,etapos,xiweight,etaweight] = GetGPV(ngaus,igaus,jgaus,xipos,etapos,xiweight,etaweight)
%% get sampling point location and associated weight factor

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 17 2017
%  Computational Structural Analysis Spring 2017

global posgp weigp

if ngaus==2
  xipos  = posgp((igaus-1)*ngaus+jgaus,1);
  etapos = posgp((igaus-1)*ngaus+jgaus,2);
  xiweight  = weigp((igaus-1)*ngaus+jgaus,1);
  etaweight = weigp((igaus-1)*ngaus+jgaus,2);
end
