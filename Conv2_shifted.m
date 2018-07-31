function [P,CC]=Conv2_shifted(t,hw,k1,fw,gw,bin)
%computes the convolution of the vectors fw and gw at the points in t.  
%In case k1 is nonzero, k1 initial values in fw are not passed. 
% To compenstate for this, C(i) approximates the colvolution of fw and gw 
% at x(i+k1), where x=0:hw:max(t).  
%Note that the convolution evaluated at x<x(k1) is zero.    

n=length(t);

%Compute fw*g %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
C=conv(fw,gw)*hw;
N=length(gw);
% only the first N-k1 elements of the convolution are valid
C=C(1:N-k1);
CC=zeros(N,1);
CC(k1+1:N)=C;
I=zeros(n,1);
% P gives the probability of each observation under the windowed
% convolution.
P=realmin*ones(n,1);
%%%%%%%%%%%%%%%%%%


%number of steps to go back, in order to go back .1 time units

if strcmp(bin,'yes')
goback=.1/hw;
goback=round(goback);
for i=1:n
    %find index of element of x that is closest to t(i)
    %[~,I(i)]=min(abs(t(i)-xw));
    I(i)=round(t(i)/hw)+1;
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
            if max(I_vector)>length(C)
                I_vector
            end
        P(i)=sum(C(I_vector))*hw;
    
    end
%toc
end
end
if strcmp(bin,'no')
    for i=1:n
        %should it be +1
        I(i)=round(t(i)/hw)+1;
    %find index of element of x that is closest to t(i)
        if t(i)>0 && I(i)>k1
            if (I(i)-k1)>length(C)
            I(i);
            end
    
            P(i)=C(I(i)-k1);
        end
    end
end

end