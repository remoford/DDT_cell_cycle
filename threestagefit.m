function [pd_max,max_ld]=threestagefit(data,TolFun,TolX,bin)

%num = length(data);
C1 = mean(data);
C2 = var(data);
%C3 = sum((data-C1).^3)/(length(data));

%Fit three-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN FUNCTION FIT_THREESTAGE
    [P]=threestage_seeds(C1,C2)

        % optimize parameters
        pd=zeros(length(P),6);
        ld = NaN*ones(length(P),1);
        fminsearch_options = optimset('TolFun',TolFun, 'TolX', TolX,'MaxFunEvals',10000,'Display','final');
        %options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        for i=1:length(P)
            x0=P(i,:);
            f=@(x,m1,s1,m2,s2,m3,s3)convolv_3invG_Dirac_option(x,m1,s1,m2,s2,m3,s3,bin);
            myll=@(params)loglikelihood(data, f, 6, params);
            objfun=@(params)penalize(myll, 6, params, [realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax]);
            [p,l]=fminsearch(objfun,x0,fminsearch_options)
            pd(i,:)=p;
            ld(i)=l;
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=min(ld);
        max_ld=-max_ld;
        pd_max = pd(row_ld,:);
        
    % END FUNCTION FIT_THREESTAGE

end