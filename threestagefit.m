function [pd_max,max_ld]=threestagefit(data)

%num = length(data);
C1 = mean(data);
C2 = var(data);
%C3 = sum((data-C1).^3)/(length(data));

%Fit three-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN FUNCTION FIT_THREESTAGE
    
        % prepare statistical parameters
        vry = [0.1 0.7]';
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        
        % prepare parameter seeds
        pcomb = allcomb(m,s);
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        id = allcomb(1:length(pcomb),1:length(pcomb),1:length(pcomb));
        id = sort(id,2);
        id = unique(id,'rows');
        P = zeros(length(id),6);
        for ii = 1:length(id)
            P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)},pcell{id(ii,3)}];
        end

        % optimize parameters
        pd=zeros(length(P),6);
        ld = NaN*ones(length(P),1);
        fminsearch_options = optimset('TolFun',TolFun, 'TolX', TolX);
        %options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        for i=1:length(P)
            x0=P(i,:);
            f=@(x,m1,s1,m2,s2,m3,s3)convolv_3invG_adapt_window(t,m1,s1,m2,s2,m3,s3);
            myll=@(params)loglikelihood(data, f, 4, params);
            objfun=@(params)penalize(myll, 6, params, [realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax]);
            p=fminsearch(objfun,x0,fminsearch_options);
            pd(i,:)=p;
            [l,hp(i)]=convolv_3invG_adapt_window(t,m1,s1,m2,s2,m3,s3);
            l=sum(log(l));
            ld(i)=l
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld_true);
        pd_max = pd(row_ld,:);
        
    % END FUNCTION FIT_THREESTAGE

end