function R=bisect(f,a,b,value)
%This code finds a point at which f is greater than value;
%x0 is a vector containing two initial points at which f is evaluated.
%it is assumed that f(a)<value and f(b)>value;
%the algorithm estimates the first point on the segment from x0(1) to x0(2) at
%which f>value.
eps_x=1;

if f(a)>value
    R=NaN;
else

while eps_x>10^(-10)
    c=(b+a)/2;
 
    if f(c)>value
    b=c;
    
    end
    if f(c)<=value
    a=c;
    
    end
    eps_x=abs(b-a)/b;
end

R=b;
end
    
end

