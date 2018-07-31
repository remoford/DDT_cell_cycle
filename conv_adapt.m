function [P,C,E,h]=conv_adapt(t,f1,f2,LL1,LL2,RL1,RL2,h0,EType,bin)
%convolves the functions f1 and f2 with support detrmined by LL and RL
%adapts the mesh size to achieve given error type

%when t is a regular grid, h0 should be the grid size.

%EB=bound on the error

E=0;

if strcmp(EType,'relLL')
%Determines bound of error in relative LL, smaller = less error.
eps=.5;
EB=min(-log(1-eps),log(1+eps));
end

if strcmp(EType,'abspdf')
EB=10^(-2);
end

if strcmp(EType,'none')
EB=0;
end

if strcmp(EType,'integral')
EB=10^(-6);
end

xw=LL2:h0:round(RL2/h0)*h0;

%gives the index of the last component of xw that is outside the support of f2.

k1=find(xw<LL1,1,'last');
w=xw(k1+1):h0:round(RL1/h0)*h0;

%begin function conv_adapt
%inputs: f1 f1: functions to be convolved
%        LL1, RL1: left and right limits of support of f1
%       LL2, RL2: left and right limits of support of f2
%       h0=initial mess for convolution
%       EType: error type
%       We assume the LL1>=LL2 and the answer is normalized to LL2.
%

f1vector=f1(w);

f2vector=f2(xw);

if length(f1vector)>length(f2vector)
        length(f1vector)
end


    
[P,C]=Conv2_shifted(t,h0,k1,f1vector,f2vector,bin);

%Only adapt if there is a bound on the error
if strcmp(EType,'none')==0
E=EB+1;
end

h=h0;
    while E>EB
    
    h=h*.5;
    
    xw=LL2:h:round(RL2/h)*h;
    k1=find(xw<LL1,1,'last');
    w=xw(k1+1):h:round(RL1/h)*h;
    
    f1vector=f1(w);

    f2vector=f2(xw);
    
    P0=P;
    
    [P,C]=Conv2_shifted(t,h,k1,f1vector,f2vector,bin);
 
    [E,problem_P]=ComputeError(P0,P,EType,EB,h0);
    
    if h<10^(-4)
        h
    end
    
    end