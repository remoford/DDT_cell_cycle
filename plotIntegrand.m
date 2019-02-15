function Q = plotIntegrand(b,d,scale)

a=1;              % parameters for our distributions
%b=1/50;
c=1;
%d=1/50;

% evaluation options
tBinLeft=1.25;      % binning boundries
tBinRight=1.5;
h=0.01;             % underlying mesh size
H=h*(2^scale)*10               % mesh size for riemann integration
lim=3;              % mesh boundry
threshold=0.1;      % value of window threshold
BW=0.2;             % bin width

% display options
maxEnable=1;        % show the maximum point
binEnable=0;        % show bin boundries
threshEnable=0;     % show threshold floor
countourEnable=0;   % show derivative contours (look underneath!)
inputEnable=0;      % show the input functions
integrandEnable=0;  % show the integrand surface
riemannEnable=1;    % show the riemann interpolated surface
errorEnable=0;      % show the error when using riemann
tauDerivEnable=0;   % show the derivative in the tau dimension
tDerivEnable=0;     % show the derivative in the t dimension
gradEnable=0;
pdfEnable=0;
pdfRiemannEnable=0;
pdfErrorEnable=0;
PEnable=0;
PRiemannEnable=0;
PErrorEnable=0;

% initialize our products
[tau,t] = meshgrid(0:h:lim, 0:h:lim);
Izeros=zeros(lim/h+1);
I=Izeros;
R=Izeros;
Iprime=Izeros;
IPrime=Izeros;
F1=Izeros;
F2=Izeros;

% loop by t
for T=0:h:lim
   tIdx = round(T/h)+1;
   
   % get the value of the integrand
   [myIntegrand, f1, f2] = integrand(0:h:lim,a,b,c,d,T);
   I(tIdx,:) = myIntegrand';
   
   % take the derivative in tau
   Iprime(tIdx,:) = gradient(myIntegrand');
   
   % save the input distributions
   F1(tIdx,:) = f1';
   F2(tIdx,:) = f2';
end

% loop by tau
for Tau=0:h:lim
    tauIdx = round(Tau/h)+1;
    
    % take the derivative in t
    IPrime(:, tauIdx) = gradient(I(:,tauIdx));
end

% loop pointwise
for T=0:h:lim
    for Tau=0:h:lim
        tIdx = round(T/h)+1;
        tauIdx = round(Tau/h)+1;
       
        % calculate indicies of riemann sum's interpolation point
        tRIdx=round(tIdx/H)*H+1;
        tauRIdx=round(tauIdx/H)*H+1;
        
        % set the value to the interpolated point's
        R(tIdx,tauIdx) = I(tRIdx,tauRIdx);
    end
end

% the error is the difference between the riemann interpolant and the true
% function
Error=R-I;

% loop by t
for T=0:h:lim
    tIdx = round(T/h)+1;
    pdf(tIdx,:) = trapz(I(tIdx,:));
    pdfRiemann(tIdx,:) = trapz(R(tIdx,:));
    pdfError(tIdx,:) = trapz(Error(tIdx,:));
    
end
Q = trapz(pdf)/(H*H);



%find the maximum point in the integrand
IMax = max(I(:));
[tauMax,tMax]=find(I==IMax);


ondemand = @(tau, k) integrand(tau, a, b, c, d, k);

[gradx, grady] = gradient(I);

%P=zeros(lim/BW,1);
%PRiemann=zeros(lim/BW,1);
%PError=zeros(lim/BW,1);
for B=0:BW:lim
   line([B B], [0 max(pdf)]);

   bIdx = round(B/h)+1;
   bIdxPrev = round(((bIdx-1)*h-BW)/h)+1;
   
   if bIdxPrev < 1
       bIdxPrev=1;
   end

   if B~=0
       pdfRange = pdf(bIdxPrev:bIdx);
       P(round(B/BW)) = trapz(pdfRange);
       
       pdfRiemannRange = pdfRiemann(bIdxPrev:bIdx);
       PRiemann(round(B/BW)) = trapz(pdfRiemannRange);
       
       pdfErrorRange = pdfError(bIdxPrev:bIdx);
       PError(round(B/BW)) = trapz(pdfErrorRange);
   end
end 

clf
hold on;

%draw the maximum point
if maxEnable==1
scatter3(tMax*h,tauMax*h,IMax,'r*');
end

%draw the threshold plane - some plots are better with this off!
if threshEnable==1
thresholdPlane=ones(lim/h+1)*threshold;
mesh(tau, t, thresholdPlane,'edgecolor','r');
end

% this helps to show the boundries of a bin
if binEnable==1
y=0:h:lim;
x = (ones(lim/h+1,1).*tBinLeft)';
z=I(tBinLeft/h,:);
plot3(y, x, z, 'LineWidth', 4)
x = (ones(lim/h+1,1).*tBinRight)';
z=I(tBinRight/h,:);
plot3(y, x, z, 'LineWidth', 4)
end

if inputEnable==1
title('Input Distributions - these are multiplied to give the integrand');
mesh(tau, t, F1);
mesh(tau, t, F2);
end

if integrandEnable==1
title('Integrand - integrate across tau from 0 to inf to get probability density then integrate across t inside the bin to get probability');
mesh(tau, t, I);
end

if riemannEnable==1
title('Left-handed Riemann sum - integrate across tau from 0 to inf to get probability density then integrate across t inside the bin to get probability');
mesh(tau, t, R);
end

if errorEnable==1
title('Error between Riemann sum and Integrand');
mesh(tau, t, Error);
end

if tauDerivEnable==1
title('Derivative in tau');
mesh(tau, t, Iprime);
end

if tDerivEnable==1
title('Derivative in t');
mesh(tau, t, IPrime);
end

if gradEnable==1
    title('Gradient');
    quiver(gradx,grady);
end

if ( pdfEnable==1) || (pdfRiemannEnable==1) || (pdfErrorEnable==1)
    xlabel('t');
    ylabel('Probability density');
    legend('pdf','riemann','error');
else
    xlabel('tau');
    ylabel('t');
end


if pdfEnable==1
   plot(t, pdf, 'g');   
end

if PEnable==1
   bar(P,'FaceAlpha',0,'EdgeColor','green'); 
end

if PRiemannEnable==1
   bar(PRiemann,'FaceAlpha',0,'EdgeColor','blue'); 
end

if PErrorEnable==1
   bar(PError,'FaceAlpha',0,'EdgeColor','red'); 
end

if pdfRiemannEnable==1   
   plot(t, pdfRiemann, 'b'); 
end

if pdfErrorEnable==1
   plot(t, pdfError, 'r'); 
end

if countourEnable==1
title('Basis derivative countours - note intersection is the integrand maximum');
meshc(tau, t, Iprime,0.001);
meshc(tau, t, IPrime,0.001);
end

annotation('textbox',[.05 .05, .3 .3],'String',num2str(Q),'FitBoxToText','on');

drawnow

end