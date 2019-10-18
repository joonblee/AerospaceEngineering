%
%  Newmark_applied_forces.m  ver 1.1  January 9, 2014
%
%
function[FP]=Newmark_applied_forces(t,ndof,iu)
%
nt=length(t);
dt=(t(nt)-t(1))/(nt-1);
tstart=t(1);
tend=t(nt);
%
FP=zeros(nt,ndof);
%
disp(' ');
out1=sprintf(' Enter the number of dofs with applied force. (maximum = %d) ',ndof);
disp(out1)
naf=input(' ');
if(naf>ndof)
    naf=ndof;
end
%
disp(' ');
if(iu==1)
    disp(' Each input file must have two columns: time(sec) & force(lbf) ');
else
    disp(' Each input file must have two columns: time(sec) & force(N) ');    
end    
%
for ijk=1:naf
%
    clear FI;
    clear yi;
    clear max;
    clear L; 
%
    disp(' ');
    out1=sprintf(' Enter force array %d:  ',ijk);
    FS = input(out1,'s');
    FI=evalin('caller',FS);
%        
%    
    disp(' Enter the dof at which this force is applied.   ')
    fdof=input(' ');
%      
%  Interpolate the force
%
    L=length(FI(:,1));
    f=zeros(nt,1);
    jmin=1;
    for i=1:nt
        f(i)=0;
        for j=jmin:L
            if(t(i)==FI(j,1))
                f(i)=FI(j,2);
                jmin=j;
                break;
            end
            if(j>=2 && FI(j-1,1) < t(i) && t(i) < FI(j,1) )
                    l=FI(j,1)-FI(j-1,1);
                    x=t(i)-FI(j-1,1);
                    c2=x/l;
                    c1=1-c2;
                    f(i)=c1*FI(j-1,2)+ c2*FI(j,2);
                    jmin=j;
                break;
            end
        end
    end
%
    FP(:,fdof)=f;
%
end