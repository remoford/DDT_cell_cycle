function R=bisect(f,a,b)
%This code finds a point at which f is greater than real min;
%x0 is a vector containing two initial points at which f is evaluated.
%it is assumed that f(a)<0 and f(b)>0;
%the algorithm estimates the first point on the segment from x0(1) to x0(2) at
%which f>realmin.
eps_x=1;

if f(a)>0
    R=NaN;
else

while eps_x>10^(-10)
    c=(b+a)/2;
 
    if f(c)>realmin
    b=c;
    
    end
    if f(c)<=realmin
    a=c;
    
    end
    eps_x=abs(b-a)/b;
end

R=b;
end
    
end

