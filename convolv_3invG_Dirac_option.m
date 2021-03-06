%this function evaluates the convolution of three inverse gaussian
%distributions at vector t
function [P,h,flag,E]=convolv_3invG_Dirac_option(t,m1,s1,m2,s2,m3,s3,bin,order)

h0=.025;
timing_output=0;
if timing_output == 1
    tic
end

flag=0;
E=Inf;
eps=.01;

m=[m1 m2 m3];
s=[s1 s2 s3];
m=max(0,m);
s=max(0,s);
v=(s.^2)./(m.^3);
[v,I]=sort(v);
sd=v.^.5;
m=m(I);
s=s(I); 

if sd(1)<0
    
    check2 = TailMass2(m,s,eps,sd);
    
    if check2<=eps/3
        
        flag=1;
        l=1/m(1);
        
    
        [P,h,flag1,E]=convolv_2invG_Dirac_option(t-l,m(2),s(2),m(3),s(3),bin,h0,'relLL');
        
        flag=flag+flag1;
    
    else
    %The next line will return the convolution of the inverse Gaussians
    %with parameters 4-7, with the inverse Gaussian with parameters 2-3.
    %We think it will be faster to first convolve the most concentrated
    %distributions and then take the convolution of that with the remaining
    %distribution.  This way the most concentrated distribution is only
    %required to meet the integral error threshold.  
        [P,h,~,E] = convolv_3invG(t,m(1),s(1),m(2),s(2),m(3),s(3),bin,h0,'relLL');
    end
        
%otherwise, perform the convolution

else
    if order==3
        [P,h,~,E] = convolv_3invG(t,m(3),s(3),m(1),s(1),m(2),s(2),bin,h0,'relLL');
    end
    if order==1
        [P,h,~,E] = convolv_3invG(t,m(1),s(1),m(2),s(2),m(3),s(3),bin,h0,'relLL');
    end
    

end
end