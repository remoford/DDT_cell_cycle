function check2 = TailMass2(m,s,eps,sd)

% BEGIN FUNCTION TailMass
        % Input parameters: m, s, eps, T2, sd,
        % Outputs: check2
        
%T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));
T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2)));
T3=(1/m(3))*((1+(9/4)*(s(3)^4/m(3)^2))^.5-(3/2)*(s(3)^2/m(3)));

           gp=min(gp_max(m(2),s(2)),gp_max(m(3),s(3)));
   
    %define the radius, r, of a small interval over which the second pdf is
    %nearly constant and close to zero for t<r. 
    T=max(T2,T3);
    r=min(eps/(3*gp),T);
    
    checkval=min(onestagepdf2(r,m(2),s(2)),onestagepdf2(r,m(3),s(3)));
    while checkval>=eps/2
        r=r/2;
        checkval=min(onestagepdf2(r,m(2),s(2)),onestagepdf2(r,m(3),s(3))) ;
    end
    
    %estimate the maximum value of the convolution of the second two pdfs
    gm=min(onestagepdf2(T2,m(2),s(2)), onestagepdf2(T3,m(3),s(3)));
    
    
    %get the average value of the first pdf.  This is the point at which
    %its mass is concentrated.
    nu=1/m(1);
    
    %get upper limit of integral for approximating
    %int_{r+nu}^{infty}f(s)ds.
    Tu=max(100,nu+r+1000*sd(1));
    
        
        checkerror=100;
        hh=.001;
        check1=.001*sum(onestagepdf2((0:.001:nu-r),m(1),s(1)))+.001*sum(onestagepdf2((nu+r:.001:Tu),m(1),s(1)));
        while checkerror>10^(-4)
            hh=.5*hh;
            ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1)));
            checkerror=abs(check1-ck1);
            check1=ck1;
        end
   check2=gm*check1;