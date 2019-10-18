function [SFD] = GetSFD(SFD,xipos,etapos)
%% compute shape function derivatives

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 17 2017
%  Computational Structural Analysis Spring 2017

% xi-derivatives d_N/d_xi
SFD(1,1) = -0.25*(1-etapos);
SFD(1,2) =  0.25*(1-etapos);
SFD(1,3) =  0.25*(1+etapos);
SFD(1,4) = -0.25*(1+etapos);

% eta-derivative d_N/d_eta
SFD(2,1) = -0.25*(1-xipos);
SFD(2,2) = -0.25*(1+xipos);
SFD(2,3) =  0.25*(1+xipos);
SFD(2,4) =  0.25*(1-xipos);

return
