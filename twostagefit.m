function [pd_max,max_ld]=twostagefit(data,TolFun,TolX)

%num = length(data);
C1 = mean(data);
C2 = var(data);
%C3 = sum((data-C1).^3)/(length(data));


%Fit two-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % BEGIN FUNCTION FIT_TWOSTAGE
    
%         % prepare statistical parameters
%         vry = [.25 .5 .75]';
%         % vrys=[.01 1 10]';
%         % [x0, l_moments, min_res] = moments_method_2stage(data);
%         % m1 = x0(1)*vry; 
%         % s1 = x0(2)*vrys;
%         % m2 = x0(3)*vry;
%         % s2 = x0(4)*vrys;
%         c1 = C1*vry;
%         c2 = C2*vry;
%         m = 1./c1;
%         s = (c2./c1.^3).^0.5;
%         N = length(vry);
%         
% 
%         % prepare parameter seeds
%         
%         %get all pairs of the form [m(i),s(j)]
%         %these pairs represent all possible unique 
%         %parameter choices for each part of the cell
%         %cycle.  
%         pcomb = allcomb(m,s);
%         %place paramter pairs into a cell.  The parameters choices for each part
%         %are now indexed
%         pcell = cell(length(pcomb),1);
%         for i = 1:length(pcomb)
%             pcell{i} = pcomb(i,:);
%         end
%         %get all pairs of indices for the parameter 
%         %choices for each part of the cycle to get all 
%         %parameter choices for the entire cycle
%         id = allcomb(1:length(pcomb),1:length(pcomb));
%         %sort the pairs in ascending order.  
%         %This equates choices of the form [i,j] and [j,i].
%         id = sort(id,2);
%         %remove repeats
%         id = unique(id,'rows');
%         %create a matrix of unique parameter choices for the cell cycle
%         P = zeros(length(id),4);
%         for ii = 1:length(id)
%             P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}];
%         end

        [P]=twostage_seeds(C1,C2)
        % optimize parameters
        pd=zeros(length(P(:,1)),4);
        ld = NaN*ones(length(P),1);
        %options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        fminsearch_options = optimset('TolFun',TolFun, 'TolX', TolX);
        flag=zeros(length(P),1);
        %confint=zeros(2,4,length(id));
        for i=1:length(P)  
            x0 = P(i,:);
            %f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_nov(x,m1,s1,m2,s2,.01);
            f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_window(x,m1,s1,m2,s2);
            myll=@(params)loglikelihood(data, f, 4, params);
            objfun=@(params)penalize(myll, 4, params, [realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax]);
            p=fminsearch(objfun,x0,fminsearch_options);
        
%             [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options)
            pd(i,:)=p
%             confint(:,:,i)=conf1;
            [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_window(data,p(1),p(2),p(3),p(4));
            l=sum(log(l));
            ld(i)=l    

        end
        
       
        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld);
        pd_max = pd(row_ld,:);
        %confint_max=confint(:,:,row_ld);
    % END FUNCTION FIT_TWOSTAGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end