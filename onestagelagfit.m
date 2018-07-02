function [pd_max,max_ld]=onestagelagfit(data)

num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));
%Fit one-stage model with lag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN FUNCTION FIT_ONESTAGELAG
    
        % prepare statistical parameters
        mu = C3/(3*C2^2);
        sigma = (C3^3/(27*C2^5))^.5;
        lag = C1-3*C2^2/C3;
        vryv=[0.5 1 2];
        vrym=[.25 .5 .75];
        m = mu*vrym;
        s = sigma*vryv;
        lag = lag*vrym;
        N = length(vrym);
        
        % prepare parameter seeds
        pp = cell(N^3);
        for i = 1:N
            for j = 1:N
                for k = 1:N
                pp{N^2*(i-1)+N*(j-1)+k} = [m(i), s(j), lag(k)];
                end
            end
        end
        P = zeros(N^3,3);
        for ii = 1:N^3
            P(ii,:) = pp{ii};
        end
        
        % optimize parameters
        pd = zeros(N^3,3);
        ld = -realmax*ones(N^3,1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        for i=1:length(P)
            x0=P(i,:);
            [p,conf1]=mle(data,'pdf',@onestagepdf_lag,'start',x0, 'upperbound', [Inf Inf Inf],'lowerbound',[0 0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            l=onestagepdf_lag(data,p(1),p(2),p(3));
            ld(i)=sum(log(l));
        end
        
        % common to each fit, consider factoring out
        [max_ld,ind_ld]=max(ld);
        pd_max=pd(ind_ld,:);
        confint_max=confint(:,:,ind_ld);
    % END FUNCTION FIT_ONESTAGELAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end