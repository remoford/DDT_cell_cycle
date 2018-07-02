function [pd_max,max_ld]=onestagefit(data)

num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));


%Fit the one-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN FUNCTION FIT_ONESTAGE
    
        % prepare statistical variables
        mu=1/C1;
        sigma=(C2/C1^3)^.5;
        vry=.5:.5:2;
        m=mu*vry;
        s=sigma*vry;
        N=length(vry);
        
        % maybe these should be moved down into optimize parameters
        % section?
        pd=zeros(N^2,2);
        ld=-realmax*ones(N^2,1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        
        % prepare parameter seeds
        pp = cell(N^2);
        for i = 1:N
            for j = 1:N
                pp{(i-1)*N+j} = [m(i), s(j)];
            end
        end
        P = zeros(N^2,2);
        for ii = 1:N^2
            P(ii,:) = pp{ii};
        end
        
        
        % optimize parameters
        for i=1:N^2
            x0 = P(i,:);
            %[p,conf1]=mle(data,'pdf',@onestagepdf2,'start',x0, 'upperbound', [Inf Inf],'lowerbound',[0 0],'options',options)
            [p,conf1]=mle(data,'pdf',@onestagepdf_binned_adapt,'start',x0, 'upperbound', [Inf Inf],'lowerbound',[0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1;
            l=onestagepdf_binned_adapt(data,p(1),p(2));
            ld(i)=sum(log(l));
        end
        
        % common to each fit, consider factoring out
        [max_ld,ind_ld]=max(ld);
        pd_max=pd(ind_ld,:);
        confint_max=confint(:,:,ind_ld);
    % END FUNCTION FIT_ONESTAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end