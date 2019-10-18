function [xipos,etapos,xiweight,etaweight] = GetGPV(ngaus,igaus,jgaus,xipos,etapos,xiweight,etaweight)
%% get sampling point location and associated weight factor

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

global posgp weigp

%two points rule
if ngaus==2

    if igaus==1
        xipos = posgp(jgaus,1);
        etapos = posgp(jgaus,2);
        xiweight = weigp(jgaus,1);
        etaweight = weigp(jgaus,2);
    else
        xipos = posgp(igaus+jgaus,1);
        etapos = posgp(igaus+jgaus,2);
        xiweight = weigp(igaus+jgaus,1);
        etaweight = weigp(igaus+jgaus,2);
    end
    
end

% you can try other gauss quadrature rules