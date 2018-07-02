function [pd_max,max_ld]=threestagefit(data)

num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));

%Fit three-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN FUNCTION FIT_THREESTAGE
    
        % prepare statistical parameters
        vry = [0.1 0.7]';
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        N = length(vry);
        
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
        flag=zeros(length(P),1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        for i=1:length(P)
            x0=P(i,:);
            g=@(x,m1,s1,m2,s2,m3,s3)convolv_3invG_nov(x,m1,s1,m2,s2,m3,s3,.01);
            [p,conf1]=mle(data,'pdf',g,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            [l,hp(i),flag(i),E(i)]=convolv_3invG_nov(data,p(1),p(2),p(3),p(4),p(5),p(6),.01);
            l=sum(log(l));
            ld(i)=l
        end

        % we previously optimized with a larger step size, recalculate with
        % a smaller stepsize after the fact
        for i=1:length(ld)
            [l,hp_true(i),flag_true(i),E_true(i)]=convolv_3invG_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),pd(i,5),pd(i,6),.001);
            ld_true(i)=sum(log(l));
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld_true);
        pd_max = pd(row_ld,:);
        confint_max=confint(:,:,row_ld);
    % END FUNCTION FIT_THREESTAGE

end