function [P]=threestage_seeds(C1,C2,vry,twostagearg)

%C1 is the mean of the data
%C2 is the variance of the data

%This code finds parameter seeds for a threestage model that match the mean 
% and variance of the data.  
%These seeds correspond to parameter combinations wherein a single part 
%of the model accounts for 10% or 50% of the variance and 10% or 50% of the
%mean.  The remaining parts of the model account for 90% or 50% of these
%moments, respectively.  
        twoseeds=twostage_seeds(C1,C2,twostagearg);
    
        % proportions of moments represented by the first part
        if strcmp(vry,'fine')
            vry1 = [.1 .2 .3 .4 .5];
        end
        if strcmp(vry,'coarse')
            vry1 = [.25 .75];
            %vry1 = [.1 .75];
        end
            
        % proportions of moments represented by the remaining parts
        portion = abs(1-vry1);
        
        c1 = C1*vry1;
        c2 = C2*vry1;
        
        c1_comp=C1*portion;
        c2_comp=C2*portion;
        

        N1 = length(vry1);
        num2seeds=size(twoseeds,1);
        
        P=zeros((N1^2)*(num2seeds),6);
        
        
        for i=1:N1
            for j=1:N1
                [P2]=twostage_seeds(c1_comp(i),c2_comp(j),twostagearg);
                for k=1:num2seeds
                    P((N1*num2seeds)*(i-1)+(num2seeds)*(j-1)+k,:)=[1/c1(i),(c2(j)/c1(i)^3).^0.5,P2(k,1),P2(k,2),P2(k,3),P2(k,4)];
                end
            end
        end
            
end
