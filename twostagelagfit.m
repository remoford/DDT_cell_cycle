function [pd_max,max_ld]=twostagelagfit(data)

num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));

%Fit two-stage model with lag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN FUNCTION FIT_TWOSTAGELAG
    
        % prepare statistical parameters
        vry = [.2 .7]';
        % vrys=[.01 1 10]';
        % [x0, l_moments, min_res] = moments_method_2stage(data);
        % m1 = x0(1)*vry; 
        % s1 = x0(2)*vrys;
        % m2 = x0(3)*vry;
        % s2 = x0(4)*vrys;
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        l = min(data)*vry;
        N = length(vry);

        % prepare parameter seeds
        
        %get all pairs of the form [m(i),s(j)]
        %these pairs represent all possible unique 
        %parameter choices for each stochastc part of the cell
        %cycle.  
        pcomb = allcomb(m,s);
        %place paramter pairs into a cell.  The parameters choices for each part
        %are now indexed
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        %get all triples of indices for the parameter 
        %choices for each part of the cycle to get all 
        %parameter choices for the entire cycle
        id = allcomb(1:length(pcomb),1:length(pcomb));
        %sort the pairs in ascending order.  
        %This equates choices of the form [i,j] and [j,i].
        id = sort(id,2);
        %remove repeats
        id = unique(id,'rows');
        %create a matrix of unique parameter choices for the cell cycle
        P = zeros(length(id)*N,5);
        for ii = 1:length(id)
            for jj=1:N
            P((ii-1)*N+jj,:) = [pcell{id(ii,1)},pcell{id(ii,2)},l(jj)];
            end
        end

        % optimize parameters
        pd=zeros(length(P),5);
        ld = NaN*ones(length(P),1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        % if flag(i)=1 the pdf is approximated as a one-part stocahstic process
        % with a dterministic lag.
        flag=zeros(length(P),1);
        for i=1:length(P)  
            x0 = P(i,:);
             f=@(t,m1,s1,m2,s2,l)twostagepdf_lag(t,m1,s1,m2,s2,l,.01,10^(-6));
            [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0],'options',options)
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            [l,flag(i)]=twostagepdf_lag(data,p(1),p(2),p(3),p(4),p(5),.01,10^(-6));
            l=sum(log(l));
            ld(i)=l
            %check error in Riemann sum
            if flag(i)==0
                [l]=twostagepdf_lag(data,p(1),p(2),p(3),p(4), p(5),.001,10^(-6));
                l=sum(log(l));
                ld_check(i)=l
            end
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld);
        pd_max = pd(row_ld,:);
        confint_max=confint(:,:,row_ld);
    % END FUNCTION FIT_TWOSTAGELAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
