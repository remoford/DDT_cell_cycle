function [P,hw,C,E] = convolv_2invG(t,m1,s1,m2,s2,bin,h0,EType)
%This function convolves a highly concentrated inverse Gaussian (f1) with 
% a less concentrated inverse Gaussian (f2).  
% It breaks the concentrated pdf as a sum of two pdfs: f=fw+fo, 
%where the domain of fw is the support of f, and fo==0.  So that
%fw*g=f*g is evaluated using an adaptive method that
%reduced the grid size, h, used to compute the convolution to ensure the
%numerical error is small. 


%m1 and s1 correspond to the more concentrated distribution

%find the left and right limits of the window over which the pdf is
%concentrated.
LL2=0;
RL2=max(t);
[LL1,RL1]=window(m1,s1,max(t));

if LL1>=max(t)
    P=realmin*ones(1,length(data));
else       

%grid size for convolving against f over the window;
hw=(RL1-LL1)/10;
hw=min(hw,h0);

% xw=linspace(0,max(t),round(max(t)/hw)+1);
% hw=xw(2)-xw(1);

%for convenience, we would like the grid size to be of the form .1/2^k.
%This ensures the data lies on a grid point.
kk=log2(.1/hw);
kk=round(kk);
hw=.1/(2^kk);

f1=@(x)onestagepdf2(x,m1,s1);
f2=@(x)onestagepdf2(x,m2,s2);

[P,C,E]=conv_adapt(t,f1,f2,LL1,LL2,RL1,RL2,hw,EType,bin);
end

