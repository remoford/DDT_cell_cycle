function [P,hw,C,E] = convolv_3invG(t,m1,s1,m2,s2,m3,s3,bin,h0,EType)
%This function convolves a highly concentrated inverse Gaussian (f1) with 
% the convolution of two less concentrated inverse Gaussians (f2 and f3).  
% It breaks the concentrated pdf as a sum of two pdfs: f=fw+fo, 
%where the domain of fw is the support of f, and fo==0.  So that
%fw*g=f*g is evaluated using an adaptive method that
%reduced the grid size, h, used to compute the convolution to ensure the
%numerical error is small. 

%The left limit of support of f2*f3 (LL23) is >= the maximum of the left
%limits of the supports of f2 and f3.  LL23>=max(LL2, LL3).
[LL2,~]=window(m2,s2,max(t));

if LL2>=max(t)
    P=realmin*ones(1,length(data));
else 
    [LL3,~]=window(m3,s3,max(t));

if LL3>=max(t)
    P=realmin*ones(1,length(data));
else 
LL23=max(LL2,LL3);
RL23=max(t);
%m1 and s1 correspond to the more concentrated distribution

%find the left and right limits of the window over which the pdf is
%concentrated.

[LL1,RL1]=window(m1,s1,max(t));

if LL1>=max(t)
    P=realmin*ones(1,length(data));
else       

%grid size for convolving against f over the window;
hw1=(RL1-LL1)/10;
hw23=(RL23-LL23)/10;
hw=min(hw1,min(hw23,h0));

%for convenience, we would like the grid size to be of the form .1/2^k.
kk=log2(.1/hw);
kk=round(kk);
hw=.1/(2^kk);

f1=@(x)onestagepdf2(x,m1,s1);

f2f3=@(x)convolv_2invG(x,m2,s2,m3,s3,'no',hw,'integral');

[P,C,E]=conv_adapt(t,f1,f2f3,LL1,LL23,RL1,RL23,hw,EType,bin);
end
end
end
end