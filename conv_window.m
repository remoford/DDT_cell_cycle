function [P] = conv_window_old(t,m1,s1,m2,s2)
%This function convolves a highly concentrated pdf (f) with a less concentrated
%pdf (g).  It breaks the concentrated pdf as a sum of two pdfs: f=fw+fo, 
%where fw=f inside the window but is zero elsewhere and fo=f outside the
%window but is zero eleswhere.  The convolution f*g=fo*g+fw*g.  The
%convolutions fo*g and fw*g are evaluated using an adaptive method that
%reduced the grid size, h, used to compute the convolution to ensure the
%numerical error is small. 


%m1 and s1 correspond to the more concentrated distribution
h=.1;
n=length(t);
%v=(s1.^2)./(m1.^3);
%get standard deviation of concentrated part for determining the window of
%the convolution.
%sd=v^.5;
%all times at which to compute the functions and their convolution
x=0:h:max(t);
x=x';
x_fine=0:.01:max(t);
%SS is size of window.  Needs to be determined.
%SS=pickthewindow(m1,s1);
%full vector of values of f.
ff=onestagepdf2(x_fine,m1,s1);
SSL=find(ff>realmin,1,'first');
SSR=find(ff>realmin,1,'last');
SSL=x_fine(SSL);
SSR=x_fine(SSR);
%grid size for convolving against f over the window;
hw=(SSR-SSL);
hw=min(hw,.1);
%vector of values at which to compute g for convolving against fw
xw=0:hw:max(t);
%location of mode of concentrated distribution
%T1=(1/m1)*((1+(9/4)*(s1^4/m1^2))^.5-(3/2)*(s1^2/m1));
%window over which the first pdf is highly concentrated
w=SSL:hw:SSR;
%vector of values of highly concentrated pdf over window
fw=onestagepdf2(w,m1,s1);
%indicates if an x value is outside of the window
%w_indicator=@(x)((x<SSL) + (x>SSR));
%gives vector of vales of highly concentrated pdf outside of window.  
%fo_fun=@(x)onestagepdf2(x,m1,s1).*w_indicator(x);
%vector of values of less concentrated pdf everywhere.
%g=onestagepdf2(x,m2,s2);
gw=onestagepdf2(xw,m2,s2);
%vector of vales of highly concentrated pdf outside of window. 
fo=fo_fun(x);

%Compute fo*g %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=conv(g,fo)*h;
N=length(g);
% only the first N elements of the convolution are valid
C=C(1:N);
%C(i) approximates f*g(x(i)) as a left-hand Riemann sum.
I=zeros(n,1);
P=zeros(n,1);
for i=1:n
    %find element of x that equals t(i).  Note with a step size of .01 all
    %of the data will fall on a grid point.
    [~,I(i)]=min(abs(t(i)-x));
    %If t(i)<0 the probability is set to zero, otherwise the
    %probability is approxiated as a value from the vector x.
    
    %number of steps to go back in order to move back .1 time units
    %because .1 is the size of the data bin FIX ME!!!!!!!
    goback=.1/h;
    if t(i)>0 && I(i)>=10
        %because the data has limited resolution
        %the probability is the integral of the pdf over the possible true 
        %IMT values.  We approximate this as a right hand Riemann sum.
        I_vector=(I(i)-goback+1):1:I(i);
        P(i)=sum(C(I_vector))*h;
    else
        P(i)=realmin;
    end
end
%toc
P0=max(realmin,P);
logP0=sum(log(P0));
%Set the initial error to be large so that the numerical convolution will
%be evaluated with at least two step sizes.
E=Inf;
%EB=bound on the error
EB=min(-log(1-.2),log(1+.2));
h1=.5*h;
while E>=EB
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
    %because we have halved the step size, the index of an observation is
    %double what it was previously.
    I=2*I-1;
    goback=goback*2;
    P=zeros(n,1);
    for i=1:n
        if t(i)>0 && I(i)>goback
            I_vector=(I(i)-goback+1):1:I(i);
            P(i)=sum(C(I_vector))*h1;
            
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
    h1=.5*h1;

end
Po=P0;
%Compute fw*g %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
C=conv(fw,gw)*hw;
N=length(gw);
% only the first N elements of the convolution are valid
C=C(1:N);
I=zeros(n,1);
% P gives the probability of each observation under the wondowed
% convolution.
P=zeros(n,1);
%gives the index of the last component of xw that is outside the window. 
k1=find(xw<SSL,1,'last');
%number of steps to go back, in order to go back .1 time units
goback=.1/xw(2);
for i=1:n
    %find element of x that is closest to t(i)
    [~,I(i)]=min(abs(t(i)-xw));
    %If t(i)<=0 the probability of an observation is set to 0.
    %Also if the index of an IMT (I(i)) is less than the first index in the 
    %window, the probability of an observation is 0.
    %Note that if t(i)<=0 means the index of the observation is 1, so such 
    %an IMT is also excluded for being outside the window. This means the
    %first condition is probably redundant.  
    if t(i)>0 && I(i)>k1
        I_vector=(I(i)-goback+1):1:I(i);
        I_vector=I_vector-k1;
        
        if 0 >= min(I_vector)
            %fprintf("OH NOES!!!!!!!\n")
            
            I_vector=I_vector(I_vector>0);
            
        end
        P(i)=sum(C(I_vector))*hw;
    else
        P(i)=0;
    end
end
%toc
%for computing the error, we just consider those observations that have
%nonzero probability when the convolution is against the windowed f.
%Note that the observed IMT must be greater than k1, for the probability 
%to be nonzero when the colvolution is against the fw.  
P0=P(I(i)>k1);
logP0=sum(log(P0));
%Set the initial error to be large so that the numerical convolution will
%be evaluated with at least two step sizes.
E=Inf;
%EB=bound on the error
EB=min(-log(1-.2),log(1+.2));
h1=.5*hw;
while E>=EB
    goback=2*goback;
    xw=0:h1:max(t);
    xw=xw';
    %window over which the first pdf is highly concentrated
    w=SSL:h1:SSR;
    %vector of values of highly concentrated pdf over window
    fw=onestagepdf2(w,m1,s1);
    %vector of values of less concentrated pdf everywhere.
    gw=onestagepdf2(xw,m2,s2);
    % BEGIN FUNCTION DOTHECONVOLUTION_INNER
    % Input parameters: z, y, h1, n, t, i, I, x
    % Outputs: logP1
    % find the discrete convolution of the vectors y and 
    % the (i-1)th element of v approximates the convolution of the pdfs 
    % over [0, x(i)] as a left-hand Riemann sum.
    C=conv(fw,gw)*h1;
    N=length(gw);
    % only the first N elements of the convolution are valid
    C=C(1:N);
    I=2*I-1;
    P=zeros(n,1);
    %first and last indices that give x values in the window
    k1=find(x<SSL,1,'last');
   for i=1:n
    %find element of x that is closest to t(i)
    [~,I(i)]=min(abs(t(i)-xw));
    %If t(i)<0 the probability is set to zero, otherwise the
    %probability is approxiated as a value from the vector x.
    if t(i)>0 && I(i)>k1
        
        %fprintf("goback: %f\n",goback);
        
        I_vector=(I(i)-goback+1):1:I(i);
        I_vector=I_vector-k1;
        
        if 0 >= min(I_vector)
            %fprintf("OH NOES!!!!!!!\n")
            
            I_vector=I_vector(I_vector>0);
            
        end
        
        %fprintf("I_vector-k1:\n");
        %I_vector-k1
        
        
        P(i)=sum(C(I_vector))*h1;
    else
        P(i)=0;
    end
    end
    %toc
    P1=P(I(i)>k1);
    logP1=sum(log(P1));
    % END FUNCTION DOTHECONVOLUTION_INNER
    E=abs(logP1-logP0);
    logP0=logP1;

    h1=.5*h1;
  
    
end
Pw=P;
P=Pw+Po;
end

