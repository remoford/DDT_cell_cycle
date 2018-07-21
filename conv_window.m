function [P] = conv_window(t,m1,s1,m2,s2,bin,h0)
%This function convolves a highly concentrated pdf (f) with a less concentrated
%pdf (g).  It breaks the concentrated pdf as a sum of two pdfs: f=fw+fo, 
%where fw=f inside the window but is zero elsewhere and fo<realmin outside the
%window and is zero eleswhere.  The convolution f*g=fo*g+fw*g.  The
%convolution fw*g is evaluated using an adaptive method that
%reduced the grid size, h, used to compute the convolution to ensure the
%numerical error is small. 

eps=.5;
%EB=bound on the error
EB=min(-log(1-eps),log(1+eps));

%m1 and s1 correspond to the more concentrated distribution

%x_fine=0:.01:max(t);
%SS is size of window.  Needs to be determined.
%SS=pickthewindow(m1,s1);
%full vector of values of f.
%ff=onestagepdf2(x_fine,m1,s1);
T1=(1/m1)*((1+(9/4)*(s1^4/m1^2))^.5-(3/2)*(s1^2/m1));
ff=@(x)(onestagepdf2(x,m1,s1));
SSL=bisect(ff,0,T1);
if SSL>=max(t)
    P=realmin*ones(1,length(data));
else
    if max(t)<T1
        SSR=max(t);
    else
    SSR=bisect(ff,max(t),T1);
    if isnan(SSR)==1
        SSR=max(t);
    end
    end

        
% SSL=find(ff>realmin,1,'first');
% SSR=find(ff>realmin,1,'last');
% SSL=x_fine(SSL);
% SSR=x_fine(SSR);
%grid size for convolving against f over the window;
hw=(SSR-SSL);
hw=min(hw,h0);
xw=0:hw:max(t);
w=SSL:hw:SSR;

fw=onestagepdf2(w,m1,s1);

gw=onestagepdf2(xw,m2,s2);

%gives the index of the last component of xw that is outside the window.
%the index of the smallest possible intermitotic time.
k1=find(xw<SSL,1,'last');
    
P=Conv2_window(t,hw,k1,fw,gw,bin);

logP=sum(log(P));



E=EB+1;

    while E>=EB
    
    hw=hw*.5
    
    xw=0:hw:max(t);
    w=SSL:hw:SSR;

    fw=onestagepdf2(w,m1,s1);

    gw=onestagepdf2(xw,m2,s2);
    
    k1=find(xw<SSL,1,'last');
    
    logP0=logP;
    
    P=Conv2_window(t,hw,k1,fw,gw,bin);
    
    logP=sum(log(P));
 
    E=abs(logP-logP0);
    
    if hw<10^(-4)
        hw
    end
    
    end
end

end

