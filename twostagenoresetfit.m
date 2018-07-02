
function [pd_max,max_ld]=twostagenoresetfit(data)

num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));

%Fit two-stage model with imperfect reset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN FUNCTION FIT_TWOSTAGE_NORESET
    
        % prepare statistical parameters
        vry = [.1 .5 .9]';
        r = [.01 .5 .99]';
        c1 = C1*vry;
        c2 = C2*vry;
        m = 1./c1;
        s = (c2./c1.^3).^0.5;
        N = length(vry);

        % prepare parameter seeds
        
        %get all pairs of the form [m(i),s(j)]
        %these pairs represent all possible unique 
        %parameter choices for each part of the cell
        %cycle.  
        pcomb = allcomb(m,s);
        %place paramter triples into a cell.  The parameter choices for each part
        %are now indexed
        pcell = cell(length(pcomb),1);
        for i = 1:length(pcomb)
            pcell{i} = pcomb(i,:);
        end
        %get all pairs of indices for the parameter 
        %choices for each part of the cycle to get all 
        %parameter choices for the entire cycle
        id = allcomb(1:length(pcomb),1:length(pcomb));
        %sort the pairs in ascending order.  
        %This equates choices of the form [i,j] and [j,i].
        id = sort(id,2);
        %remove repeats
        id = unique(id,'rows');
        %create a matrix of unique parameter choices for the cell cycle
        P = zeros(length(id)*length(r),5);
        
                
        % STUB CODE STUB CODE
        % this needs to be fixed so that there are a variety of starting
        % points for the r parameter. As is, this hard codes a starting
        % point of 0.5 which is not what we really want!
        for ii = 1:length(id)
            for jj=1:length(r)
            P(length(r)*(ii-1)+jj,:) = [pcell{id(ii,1)},pcell{id(ii,2)},r(jj)];
            end
        end

        % optimize parameters
        pd=zeros(length(P),5);
        ld = NaN*ones(length(P),1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        flag=zeros(length(id),1);
        for i=1:length(id)
            startOptimization=tic;
            x0 = P(i,:);
            fprintf("optimizing seed %d: m1=%f s1=%f m2=%f s2=%f r=%f\n", i, x0(1),x0(2),x0(3),x0(4),x0(5));
            f=@(x,m1,s1,m2,s2,r)convolv_2invG_noreset(x,m1,s1,m2,s2,r,.01);
            [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf 1],'lowerbound',[0 0 0 0 0],'options',options);
            fprintf("optimized: m1=%f s1=%f m2=%f s2=%f r=%f\n", p(1),p(2),p(3),p(4),p(5));
            pd(i,:)=p;
            confint(:,:,i)=conf1(:);
            [l,hp(i),flag(i),E(i)]=convolv_2invG_noreset(data,p(1),p(2),p(3),p(4),p(5),.01);
            l=sum(log(l));
            fprintf("log-liklihood=%f\n",l);
            if flag(i) == 1
                fprintf("used dirac delta approximation in final result\n");
            end
            ld(i)=l;
            toc(startOptimization)
            fprintf("\n");
        end
        
        % we previously optimized with a larger step size, recalculate with
        % a smaller stepsize after the fact
        ld_true=zeros(length(ld),1);
        for i=1:length(ld)
            [l,hp_true(i),flag_true(i),E_true(i)]=convolv_2invG_noreset(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),pd(i,5),.001);
            ld_true(i)=sum(log(l));
        end

        % common to each fit, consider factoring out
        [max_ld,row_ld]=max(ld_true)
        pd_max = pd(row_ld,:)
        confint_max=confint(:,:,row_ld)
    % END FUNCTION FIT_TWOSTAGE_NORESET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end