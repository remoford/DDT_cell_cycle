% This function evaluates a mixture of a convolution of two inverse Gaussian
% distributions and a single inverse Gaussian distribution at vector t.

% t is a vector of times to divide (or times to complete two parts of the
% cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1,
% s2=sigma2.

function [P,h,flag,E]=convolv_2invG_noreset(t,m1,s1,m2,s2,r,h)
% Input parameters:
% t = the list of points at which to evaluate the distribution
% m1 = mu for the first distribution
% s1 = sigma for the first distribution
% m2 = mu for the second distribution
% s2 = sigma for the second distribution
% r = probability of skipping the highly varaible cell cycle part
% h = step size

% Outputs:
% P = probability at each point corresponding to points in t
% h = final step size
% flag = indicates if we used the Dirac delta
% E = is the relative error in the likelihood of the data due to the numerical integration

% log the parameters we were called with
%fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ... 
%    m1, s1, m2, s2);


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

% find the standard deviation for each sub-distribution
sd=v.^.5;

% reorder m and s so that the sub-distribution with the smallest
% variance comes first.  So the fist part might be approximated as a Dirac delta.
[v,I]=sort(v);
m=m(I);
s=s(I);

[P2,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h);

%if the Dirac delta model is not chosen, use the inverse Gaussian for the
%probability of data, given noreset.
if flag~=1
P1=onestagepdf2(t,m(1),s(1));
end
%if the Dirac delta model is chosen, use it to compute the probability of 
% the data given no reset. 
if flag==1
P1=sum((1/m(1)-.1<=t).*(1/m(1)+.1>=t))/.2;
end

P=r*P1+(1-r)*P2;

end
