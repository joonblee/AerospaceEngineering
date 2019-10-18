function [SF] = GetSF(SF,xipos,etapos)
%% shape function

SF = 0.25 * [(1-xipos)*(1-etapos);
             (1+xipos)*(1-etapos);
             (1+xipos)*(1+etapos);
             (1-xipos)*(1+etapos)];

return
