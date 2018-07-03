function Y=onestagepdf_MATLAB(t,mu,s)
% find the value of the inverse gaussian with parameters mu and s
% at each point in t and return a list Y with the value cooresponding
% to the points in t
nu=1/mu;
lambda=1/s^2;
d = makedist('InverseGaussian',nu,lambda);

Y=pdf(d,t);

%The pdf may return values that are zero to witin machine error
%these values are also replaced by realmin
%Y=max(Y, realmin);
Y=reshape(Y,length(Y),1);

Y(isnan(Y))=0;

end