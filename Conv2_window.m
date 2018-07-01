function P=Conv2_window(SSL,SSR,hw,t,m1,s1,m2,s2)
n=length(t);
xw=0:hw:max(t);
%location of mode of concentrated distribution
%T1=(1/m1)*((1+(9/4)*(s1^4/m1^2))^.5-(3/2)*(s1^2/m1));
%window over which the first pdf is highly concentrated
w=SSL:hw:SSR;
%vector of values of highly concentrated pdf over window
fw=onestagepdf2(w,m1,s1);

gw=onestagepdf2(xw,m2,s2);


%Compute fw*g %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
C=conv(fw,gw)*hw;
N=length(gw);
% only the first N elements of the convolution are valid
C=C(1:N);
I=zeros(n,1);
% P gives the probability of each observation under the wondowed
% convolution.
P=realmin*ones(n,1);
%%%%%%%%%%%%%%%%%%
%Note that finding and integrating across the bin of each data point could
%be its own function.

%gives the index of the last component of xw that is outside the window.
%the index of the smallest possible intermitotic time.
k1=find(xw<SSL,1,'last');
%number of steps to go back, in order to go back .1 time units
goback=.1/xw(2);
for i=1:n
    %find index of element of x that is closest to t(i)
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
    
    end
%toc
end
end