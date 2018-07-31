function [P,h]=convolv_3invG_adapt_window(t,m1,s1,m2,s2,m3,s3,bin)
%onlt the first convolution is windowed
h0=.025;
Maxt=max(t);


m=[m1 m2 m3];
s=[s1 s2 s3];
m=max(0,m);
s=max(0,s);
v=(s.^2)./(m.^3);
[v,I]=sort(v);
sd=v.^.5;
m=m(I);
s=s(I);

T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));
T3=(1/m(3))*((1+(9/4)*(s(3)^4/m(3)^2))^.5-(3/2)*(s(3)^2/m(3)));


x=0:h0:Maxt;
%return the convolution of the first two over a grid at least as fine as h0;
[P,h,flag,E1,C] =convolv_2invG_Dirac_option(x,m1,s1,m2,s2,'no',.05,'abspdf');

%.1/2^kk is the grid size returned
kk=.1/h;
kk=log2(kk);
kk=round(kk);

hh=.1;

jj=.1/hh;
jj=log2(jj);
jj=round(jj);

%Maxt will be divisible by hh and h;
xhh=0:hh:Maxt;
xh=0:h:Maxt;

z=onestagepdf2(xhh,m(3),s(3));
course_index=1:2^(kk-jj):length(xh);
y=C(course_index);

eps=.5;
%EB=bound on the error
EB=min(-log(1-eps),log(1+eps));

if length(y)~=length(z)
    length(y);
end

P=Conv2_window(t,hh,0,z,y,bin);

E=EB+1;

    while E>=EB
        
        P0=P;
    
        hh=hh*.5;
    
        xhh=0:hh:Maxt;
    
        jj=.1/hh;
        jj=log2(jj);
        jj=round(jj);
    
        z=onestagepdf2(xhh,m(3),s(3));
    
        if jj<=kk
            course_index=1:2^(kk-jj):length(xh);
            y=C(course_index);
            if length(y)~=length(z)
                length(y)
            end

    
            P=Conv2_window(t,hh,0,z,y,bin);
        end
        if jj>kk
            [y]=convolv_2invG_Dirac_option(xhh,m1,s1,m2,s2,'no',2*hh,'abspdf');
            P=Conv2_window(t,hh,0,z,y,bin);
        end
 
    E=ComputeError(P0,P,'relLL',EB);
    
    end
end


