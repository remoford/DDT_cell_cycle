% This function evaluates a convolution of two inverse Gaussian
% distributions with variable initial value at vector t.  The initial value
% is uniformly distributed over an interval of radius r.


function [P,h,flag,E]=convolv_2invG_var_reset(t,m1,s1,m2,s2,r,h)
% Input parameters:
% t = the list of points at which to evaluate the distribution
% m1 = mu for the first distribution
% s1 = sigma for the first distribution
% m2 = mu for the second distribution
% s2 = sigma for the second distribution
% r = radius of interval over which reset values are uniformly distributed
% h = step size

% Outputs:
% P = probability at each point corresponding to points in t
% h = final step size
% flag = indicates if we used the Dirac delta
% E = is the relative error in the likelihood of the data due to the numerical integration

% log the parameters we were called with
fprintf('start var_reset: m1=%f s1=%f m2=%f s2=%f r=%f\n',m2,s1,m2,s2,r);

if r < 0
    fprintf("WARNING: r < 0, taking absolute value!!!\n");
    r = -r;
end
% number of points to be evaluated
n=length(t);

% store m and s for each sub-distribution in a list
m=[m1 m2];
s=[s1 s2];

% make sure we dont have negative values for m or s by replacement with
% zero if we do
m=max(m,0);
s=max(s,0);

% find the variance for both sub-distributions
v=(s.^2)./(m.^3);

% reorder m and s so that the sub-distribution with the smallest
% variance comes first.  So the fist part might be approximated as a Dirac delta.
[v,I]=sort(v);
m=m(I);
s=s(I);

%grid of possible reset values
EE=10;
hh=r;
y0=-r:hh:r;
a=1-y0;
%numerically integrate the pdg against the distribution of initial values 
PP=zeros(length(t),length(y0));
for i=1:(length(y0)-1)
[PP(:,i),~,flag,E]=convolv_2invG_adapt_nov(t,m(1),s(1),m(2)/a(i),s(2)/a(i),h);
end
%P is a vector that gives the probability of each element in t.
P=(hh/(2*r))*sum(PP,2);
logP1=sum(log(P));
 %adapt the step size in the numerical integration until the relative error 
 %in the log likelihood is small.
while EE>.001*abs(logP1)
    hh=h*.5;
    logP0=logP1;
    for i=1:(length(y0)-1)
        [PP(:,i),~,flag,E]=convolv_2invG_adapt_nov(t,m(1),s(1),m(2)/a(i),s(2)/a(i),h);
    end
    %P is a vector that gives the probability of each element in t.
    P=(hh/(2*r))*sum(PP,2);
    logP1=sum(log(P));
    E=abs(logP1-logP0);
end


%if flag == 1
%    fprintf("WARNING: Applied dirac delta approximiation, r inconsequential!\n");
%end

end
