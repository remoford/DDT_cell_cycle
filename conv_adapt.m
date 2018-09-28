function [P,C,E,h]=conv_adapt(t,f1,f2,LL1,LL2,RL1,RL2,h0,EType,bin)
%convolves the functions f1 and f2 with support detrmined by LL and RL
%adapts the mesh size to achieve given error type

%when t is a regular grid, h0 should be the grid size.

%EB=bound on the error

maxt=max(t);


if strcmp(EType,'relLL')
%Determines bound of error in relative LL, smaller = less error.
eps=.5;
EB=min(-log(1-eps),log(1+eps));
%also the integral error should be small
EB2=10^(-2);
E=EB+1;
E2=EB2+1;
end

if strcmp(EType,'abspdf')
EB=10^(-2);
EB2=0;
E=EB+1;
E2=-1;
end

if strcmp(EType,'none')
EB=0;
EB2=0;
E=-1;
E2=-1;
end

if strcmp(EType,'integral')
EB=10^(-6);
EB2=0;
E=EB+1;
E2=-1;
end

%vector spanning the support of both pdfs.
xw=0:h0:round(maxt/h0)*h0;
N=length(xw);

%gives the index of the last component of xw that is outside the support of f1.
k1=find(xw<=LL1,1,'last');
k2=find(xw<=LL2,1,'last');
w1=xw(k1):h0:round(RL1/h0)*h0;
w2=xw(k2):h0:round(RL2/h0)*h0;


f1vector=f1(w1);

f2vector=f2(w2);


    
[P,C]=Conv2_shifted(t,h0,k1+k2-1,N,f1vector,f2vector,bin);

h=h0;
    while E>EB || E2>EB2
    
    h=h*.5;
    
    xw=0:h:round(maxt/h)*h;
    N=length(xw);

    %xw(k1) is guaranteed to be just outside the support of f2, provided
    %the mesh size is >>10^(-10).

    k1=find(xw<=LL1,1,'last');
    k2=find(xw<=LL2,1,'last');
    w1=xw(k1):h:round(RL1/h)*h;
    w2=xw(k2):h:round(RL2/h)*h;


    f1vector=f1(w1);

    f2vector=f2(w2);
    
    P0=P;
    C0=C;
    
    [P,C]=Conv2_shifted(t,h,k1+k2,N,f1vector,f2vector,bin);
    %P gives the convolution at the data point, C gives the convolution on
    %the grid.
 
    [E,E2,problem_P]=ComputeError(P0,P,EType,EB,h0,C,C0);
    
    if h<10^(-4)
        h
    end
    
    end