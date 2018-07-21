function Y=onestagepdf_binned_adapt(imt,mu,s,h0)
%EB=bound on the error
h=h0;
eps=.5;
EB=min(-log(1-eps),log(1+eps));

E=EB+1;
Y0=onestagepdf_binned(imt,mu,s,h);
LogY0=sum(log(Y0));
while E>=EB
    
    h=h*.5;
    Y=onestagepdf_binned(imt,mu,s,h);
    LogY=sum(log(Y));
    E=abs(LogY-LogY0);
    LogY0=LogY;
    
end
%fprintf("gridsize=%f, E=%f, EB=%f, LL=%f\n",h,E, EB, LogY);
%fprintf("Parameters: mu=%f, sigma=%f\n\n",mu, s);
end