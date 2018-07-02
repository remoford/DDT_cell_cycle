function [pd_max,max_ld]=emgfit(data)

%Fit the EMG model 

num = length(data);
C1 = mean(data);
C2 = var(data);
C3 = sum((data-C1).^3)/(length(data));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % BEGIN FUNCTION FIT_EMG
    % Input parameters: C1, C2, data
    % Outputs: 
    
        % prepare statistical variables
        vry = [.25 .5 .75]';  
        c1=C1*vry;
        c2=C2*vry;
        %we vary the parameters so the the Gaussian and exponential parts of
        %the cell cycle are responsible for a fraction of the total mean and
        %variance in the IMT.
        lam_v=1./c1;
        mu_v=c1;
        sig_v=c2.^.5;
        N = length(vry);
        
        % prepare parameter seeds
        pp = cell(N^3);
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    pp{N^2*(i-1)+N*(j-1)+k} = [lam_v(i), sig_v(j), mu_v(k)];
                end
            end
        end
        P = zeros(N^3,3);
        for ii = 1:N^3
            P(ii,:) = pp{ii};
        end
        
        % optimize parameters
        ep = zeros(N^3,3);
        le = -realmax*ones(N^3,1);
        options = statset('MaxIter',10000, 'MaxFunEvals',10000);
        for i=1:length(P)
            x0=P(i,:);
            options = statset('MaxIter',10000, 'MaxFunEvals',10000);
            [ep(i,:),econf]=mle(data,'pdf',@emgpdf,'start',x0, 'lowerbound',[0 0 0],'upperbound',[100 100 100],'options',options);
            econfint(:,:,i)=econf(:);
            l=emgpdf(data,ep(i,1),ep(i,2),ep(i,3));
            le(i)=sum(log(l));
        end
        
        % common to each fit, consider factoring out
        [max_le,ind_le]=max(le);
        ep_max=ep(ind_le,:);
        confint_max=econfint(:,:,ind_ld);
    % END FUNCTION FIT_EMG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end