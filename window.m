function [SSL,SSR]=window(m1,s1,Maxt)
%find the left and right limits of the window over which the inverse 
%Gaussian is concentrated.  Outside this window the distribution is less
%than thresh
T1=(1/m1)*((1+(9/4)*(s1^4/m1^2))^.5-(3/2)*(s1^2/m1));
ff=@(x)(onestagepdf2(x,m1,s1));
thresh=0;
SSL=bisect(ff,0,T1,thresh);

    if Maxt<T1
        SSR=Maxt;
    else
    SSR=bisect(ff,Maxt,T1,thresh);
    if isnan(SSR)==1
        SSR=Maxt;
    end
    end
end