% This function evaluates the convolution of two inverse Gaussian
% distributions at vector t.

% If the variance in one of the distributions is very small so that the 
% distribution is close the a Dirac delta function, the convolution
% is approximated as a shifted inverse Gaussian, that is, 
% one part of the cell cycle is treated as deterministic in length.  
% In this case, the shift or lag is the average time to complete 
% the part of the cycle that has the smallest standard deviation.

% t is a vector of times to divide (or times to complete two parts of the
% cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1,
% s2=sigma2.

function [P,h,flag,E]=convolv_2invG_adapt_window(t,m1,s1,m2,s2)
h=.1;
timing_output=0;
if timing_output == 1
    tic
end

%fprintf('start 2invG: m1=%f s1=%f m2=%f s2=%f \n',m2,s1,m2,s2);
% Input parameters:
% t = the list of points at which to evaluate the distribution
% m1 = mu for the first distribution
% s1 = sigma for the first distribution
% m2 = mu for the second distribution
% s2 = sigma for the second distribution
% h = step size

% Outputs:
% P = probability at each point corresponding to points in t
% h = final step size
% flag = indicates if we used the Dirac delta
% E = is the relative error in the likelihood of the data due to the numerical integration

% log the parameters we were called with
%fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ... 
%    m1, s1, m2, s2);

% ????
flag=0;

% a constant that is used in determining if the Dirac approximation should be applied.
eps=.01;

% ????
E=Inf;

% number of points to be evaluated
n=length(t);

% store m and s for each sub-distribution in a list
m=[m1 m2];
s=[s1 s2];

% make sure we dont have negative values for m or s by replacement with
% zero if we do
m=max(m,0);
s=max(s,0);

% find the variance for both sub-distributions
v=(s.^2)./(m.^3);



% reorder m and s so that the sub-distribution with the smallest
% variance comes first.  So the fist part might be approximated as a Dirac delta.
[v,I]=sort(v);
% find the standard deviation for each sub-distribution
sd=v.^.5;
m=m(I);
s=s(I);


% T1 appears to be unused, should probably remove
% T2 is only used inside FUNCTION CHECK_APPROXIMATABLE
% and should probably be moved inside that scope
% T2 is the location of the mode of the second distribution, 
% determined by a precomputed analytic expression.
%T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));

%When called from convolv_3invG all of the components of t may be negative.
%In this case, the probability of each is set to realmin
if max(t)<=0
    P=realmin.*ones(size(t));
    P=reshape(P,length(P),1);
    
% otherwise we need to calculate P for each point in t
else
    % find the largest point in t
    Maxt=max(t);
   
    % produce a range of evenly spaced points to evaluate at between
    % 0 and Maxt with step size h. the even spacing is important when
    % calculating the convolution later
    
    %x are points at which the evaluate p2.
    x=0:h:Maxt;
    x=x';
    
    % evaluate the sub-distributions at each point in x
    y=onestagepdf2(x,m(1),s(1));
    z=onestagepdf2(x,m(2),s(2));
    
if sd(1)<=1000 && sd(1)>=.1
        
        P=conv_window(t,m(1),s(1),m(2),s(2));
    
    % if the first pdf is very concentrated, check to see if it can be
    % approximated as a point-mass distribution
end
if sd(1)<.1
    
     check2 = TailMass(m,s,eps,T2,sd);
        
        if  check2<=eps/3
            %If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1.
            %Flag is set to 1 to indicate this.
            flag=1;
            sigma=s(2);
            mu=m(2);
            %l for lag.
            l=1/m(1);
            P=onestagepdf_lag(t,mu,sigma,l);
            
            %Error in the approximation is may be too large, compute the
            %full convolution using the window.
        else
            P=conv_window(t,m(1),s(1),m(2),s(2));
        end
end
    % pdf is not very concentrated so compute the convolution directly
    if sd(1)>=1000
        
        P = Conv2(z,y,h,n,t,x);
        
        logP=sum(log(P));
       
        EB=min(-log(1-.2),log(1+.2));
        
        E=EB+1;
        
        while E>=EB
            
            h=.5*h;
            
            logP0=logP;
                
            if timing_output == 1
                fprintf('h1=%f ',h1);
                fprintf('m1=%f s1=%f m2=%f s2=%f h=%f  ',m2,s1,m2,s2,h);
            end
                x=0:h:Maxt;
                x=x';
                y=onestagepdf2(x,m(1),s(1));
                z=onestagepdf2(x,m(2),s(2));
                
                P = Conv2(z,y,h,n,t,x);
                
                logP=sum(log(P));

                E=abs(logP-logP0);
                
                if timing_output == 1
                    toc
                end
        end
           
    end
end
if timing_output == 1
    fprintf('m1=%f s1=%f m2=%f s2=%f h=%f  ',m2,s1,m2,s2,h);
    toc
    fprintf('\n\n');
end
P=reshape(P,length(P),1);
end
