function [f, f1, f2] = integrand(tau,a,b,c,d,k)




f1 = onestagepdf_MATLAB(tau,a,b);

f2 = onestagepdf_MATLAB(k-tau,c,d);

f = f1.*f2;

end