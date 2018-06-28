function Y=onestagepdf_binned_adapt(imt,mu,s)

h=.1;
E=1;
Y=onestagepdf_binned(imt,mu,s,h);

while E>.1
    h=h*.5;
    Y0=Y;
    Y=onestagepdf_binned(imt,mu,s,h);
    E=abs(Y-Y0);
end
end