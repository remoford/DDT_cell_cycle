function [P] = conv_window(t,m1,s1,m2,s2)
%This function convolves a highly concentrated pdf (f) with a less concentrated
%pdf (g).  It breaks the concentrated pdf as a sum of two pdfs: f=fw+fo, 
%where fw=f inside the window but is zero elsewhere and fo=f outside the
%window but is zero eleswhere.  The convolution f*g=fo*g+fw*g.  The
%convolutions fo*g and fw*g are evaluated using an adaptive method that
%reduced the grid size, h, used to compute the convolution to ensure the
%numerical error is small. 

%m1 and s1 correspond to the more concentrated distribution
h=.01;
n=length(t);
v=(s1.^2)./(m1.^3);
%get standard deviation of concentrated part for determining the window of
%the convolution.
sd=v^.5;
%SS is size of window.  Needs to be determined.
SS=40*sd;
%all times at which to compute the functions and their convolution
x=0:h:max(t);
x=x';
%location of mode of concentrated distribution
T1=(1/m1)*((1+(9/4)*(s1^4/m1^2))^.5-(3/2)*(s1^2/m1));
%window over which the first pdf is highly concentrated
%For now the size of the window is set to h=.01, let's pick a better size.
w=T1-SS:h:T1+SS;
%vector of values of highly concentrated pdf over window
fw=onestagepdf2(w,m1,s1);
%indicates if an x value is outside of the window
w_indicator=@(x)((x<T1-SS) + (x>T1+SS));
%gives vector of vales of highly concentrated pdf outside of window.  
fo_fun=@(x)onestagepdf2(x,m1,s1).*w_indicator(x);
%vector of values of less concentrated pdf everywhere.
g=onestagepdf2(x,m2,s2);
%vector of vales of highly concentrated pdf outside of window. 
fo=fo_fun(x);

%Compute fo*g %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=conv(g,fo)*h;
N=length(g);
% only the first N elements of the convolution are valid
C=C(1:N);
I=zeros(n,1);
P=zeros(n,1);
for i=1:n
    %find element of x that is closest to t(i)
    [~,I(i)]=min((t(i)-x).^2);
    %If t(i)<0 the probability is set to zero, otherwise the
    %probability is approxiated as a value from the vector x.
    if t(i)>0 && I(i)>1
        P(i)=C(I(i)-1);
    else
        P(i)=realmin;
    end
end
%toc
P0=max(realmin,P);
logP0=sum(log(P0));
%Set the initial error to be large so that the numerical convolution will
%be evaluated with at least two step sizes.
E=abs(logP0);
while E>=.001*abs(logP0)
    h1=.5*h;
    x=0:h1:max(t);
    x=x';
    %vector of vales of highly concentrated pdf outside of window.  
    fo=fo_fun(x);
    %vector of values of less concentrated pdf everywhere.
    g=onestagepdf2(x,m2,s2);
    % BEGIN FUNCTION DOTHECONVOLUTION_INNER
    % find the discrete convolution of the vectors g and fo 
    % the (i-1)th element of C approximates the convolution of the pdfs 
    % over [.001, x(i)] as a left-hand Riemann sum.
    C=conv(g,fo)*h1;
    N=length(g);
    % only the first N elements of the convolution are valid
    C=C(1:N);
    I=zeros(n,1);
    P=zeros(n,1);
    for i=1:n
        %find element of x that is closest to t(i)
        [~,I(i)]=min((t(i)-x).^2);
        %If t(i)<0 the probability is set to zero, otherwise the
        %probability is approximated as a value from the vector x.
        if t(i)>0 && I(i)>1
            P(i)=C(I(i)-1);
        else
            P(i)=realmin;
        end
    end
    %toc
    P1=max(realmin,P);
    logP1=sum(log(P1));
    % END FUNCTION DOTHECONVOLUTION_INNER
    E=abs(logP1-logP0);
    P0=P1;
    logP0=logP1;
    h=h1;
end
Po=P0;
%Compute fw*g %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
C=conv(fw,g)*h;
N=length(g);
% only the first N elements of the convolution are valid
C=C(1:N);
I=zeros(n,1);
% P gives the probability of each observation under the wondowed
% convolution.
P=zeros(n,1);
%gives the index of the first component of x that is inide the window. 
k1=find(x<=T1-SS,1,'last');
for i=1:n
    %find element of x that is closest to t(i)
    [~,I(i)]=min((t(i)-x).^2);
    %If t(i)<=0 the probability of an observation is set to 0.
    %Also if the index of an IMT (I(i)) is less than the first index in the 
    %window, the probability of an observation is 0.
    %Note that if t(i)<=0 means the index of the observation is 1, so such 
    %an IMT is also excluded for being outside the window. This means the
    %first condition is probably redundant.  
    if t(i)>0 && I(i)>=k1
        P(i)=C(I(i)-1-k1);
    else
        P(i)=0;
    end
end
%toc
%for computing the error, we just consider those observations that have
%nonzero probability when the convolution is against the windowed f.
%Note that the observed IMT must be greater than k1, for the probability 
%to be nonzero when the colvolution is against the fw.  
P0=P(I(i)>=k1);
logP0=sum(log(P0));
%Set the initial error to be large so that the numerical convolution will
%be evaluated with at least two step sizes.
E=abs(logP0);
while E>=.001*abs(logP0)
    h1=.5*h;
    
    x=0:h1:max(t);
    x=x';
    %window over which the first pdf is highly concentrated
    w=T1-SS:h1:T1+SS;
    %vector of values of highly concentrated pdf over window
    fw=onestagepdf2(w,m1,s1);
    %vector of values of less concentrated pdf everywhere.
    g=onestagepdf2(x,m2,s2);
    % BEGIN FUNCTION DOTHECONVOLUTION_INNER
    % Input parameters: z, y, h1, n, t, i, I, x
    % Outputs: logP1
    % find the discrete convolution of the vectors y and 
    % the (i-1)th element of v approximates the convolution of the pdfs 
    % over [0, x(i)] as a left-hand Riemann sum.
    C=conv(fw,g)*h1;
    N=length(g);
    % only the first N elements of the convolution are valid
    C=C(1:N);
    I=zeros(n,1);
    P=zeros(n,1);
    %first and last indices that give x values in the window
    k1=find(x<=T1-SS,1,'last');
   for i=1:n
    %find element of x that is closest to t(i)
    [~,I(i)]=min((t(i)-x).^2);
    %If t(i)<0 the probability is set to zero, otherwise the
    %probability is approxiated as a value from the vector x.
    if t(i)>0 && I(i)>=k1
        P(i)=C(I(i)-1-k1);
    else
        P(i)=0;
    end
    end
    %toc
    P1=P(I(i)>=k1);
    logP1=sum(log(P1));
    % END FUNCTION DOTHECONVOLUTION_INNER
    E=abs(logP1-logP0);
    logP0=logP1;
    h=h1;
    
end
Pw=P;
P=Pw+Po;
end

