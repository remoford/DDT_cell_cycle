function gpudemo(mode, numPoints)

tic;
h=10*(1/numPoints);

x=0:h:10;
if mode=='gpu'
    x=gpuArray(x);
end
x=x';

y=onestagepdf2(x,0.234,0.11);
z=onestagepdf2(x,0.5,0.2);

if mode=='gpu'
    y=gpuArray(y);
    z=gpuArray(z);
end

C=conv(z,y)*h;

figure;
hold on;
plot(y,'red');
plot(z,'blue');
plot(C,'black');

toc

end
