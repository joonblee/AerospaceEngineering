function [IDArray] = GetIDArray()
%% Get ID Array

%  Coded by Prof. Gun Jin Yun (gunjin.yun@snu.ac.kr)
%  May 15 2017
%  Computational Structural Analysis Spring 2017

global npoin ndofn 
global iffix nofix
global neq
	
% initialize of IDArray
for ipoin=1:npoin
	iposi = (ipoin-1)*ndofn;
	for idofn=1:ndofn		
		if (iffix(iposi+idofn)==0) % Active DOF
			IDArray(idofn,ipoin) = 0;            
        elseif (iffix(iposi+idofn)~=0) % Non-acticve DOF
			IDArray(idofn,ipoin) = 1;
        end
    end
end

% Generate a table of equation number
neq = 0;
for ipoin=1:npoin
	for idofn=1:ndofn
		% transfer if DOF is fixed. Otherwise, increment neq
		if (IDArray(idofn,ipoin)==0) 
			neq = neq+1;
			IDArray(idofn,ipoin) = neq;
        else 
			IDArray(idofn,ipoin) = 0;
        end
    end
end

return
