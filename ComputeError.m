function [E,E2,problem_P]=ComputeError(P0,P,EType,EB,h,C,C0)
%P0 and P are vectors that give approximation of a pdf at data points.
%If EType=relLL E is the relative error in the likelihood of the
%data
%If Etype=abspdf, E is the maximum pointwise relative error on the pdf
%values.
%If Etype=integral, P and P0 give the value of the pdf on a regular grid,
%this grid should be fine enough that the inegral across the full support
%is very close to one.  If the grid does not span the support, we only
%require the integral be less than one.
problem_P=-1;
if strcmp(EType,'relLL')
logP=sum(log(P));
logP0=sum(log(P0));
    if (logP0==-Inf && logP~=-Inf) || (logP==-Inf && logP0~=-Inf)
        E=EB+1;
    else
        if logP0==-Inf && logP==-Inf
            E=0;
        else
        E=abs(logP-logP0);
        end
    end
    %require the integral of the pdf to converge as well
    if C(length(C))==0 && max(C)~=0
    E2=abs(1-sum(C)*h);
    else
    E2=abs((sum(C)*h)-(sum(C0)*2*h));
    end
end

if strcmp(EType,'abspdf')
    %if E below satisfies E<0 the pointwise relative error in the pdf is
    %less than EB.  This formulation avoids dividing by zero.
    E=abs(P-P0)-EB*P-EB;
    [E,E_ind]=max(E);
    %Find the value of P that is responsible for the greatest error
    problem_P=P(E_ind);
    %If the bound is not met, set the error to EB+1;
    if E>0
        E=EB+1;
    end
    E2=0;
end

if strcmp(EType,'integral')
    
    if P(length(P))==0 && max(P)~=0
    E=abs(1-sum(P)*h);
    else
    E=abs((sum(P)*h)-(sum(P0)*h));
    end
    
    E2=0;
end
end


