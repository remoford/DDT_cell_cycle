function [P]=twostage_seeds(C1,C2,vry)

%C1 is the mean of the data
%C2 is the variance of the data

%This code finds parameter seeds for a twostage model that match the mean 
% and variance of the data.  
%These seeds correspond to parameter combinations wherein a single part 
%of the model accounts for 10% or 50% of the variance and 10% or 50% of the
%mean.  The remaining part of the model accounts for 90% or 50% of these
%moments, respectively.  

        % proportions of moments represented by the first part
        if strcmp(vry,'coarse')
        vry1 = [.1 .5 1.25]';
        end
        
        if strcmp(vry,'fine')
        vry1 = [.1 .2 .3 .4 .5]';
        end
    
        % proportions of moments represented by the second part
        vry2 = abs(1-vry1);
        
        c1 = C1*vry1;
        c2 = C2*vry1;
        
        c1_comp=C1*vry2;
        c2_comp=C2*vry2;
       
        N = length(vry1);
        
        P=zeros(N^2,4);
        cumulant_seeds=zeros(N^2,4);
        
        for i=1:N
            for j=1:N
            cumulant_seeds(N*(i-1)+j,:)=[c1(i),c2(j),c1_comp(i),c2_comp(j)];
            end
        end
        
        for i=1:N^2
            P(i,:)=[1/cumulant_seeds(i,1),(cumulant_seeds(i,2)/cumulant_seeds(i,1)^3).^0.5,1/cumulant_seeds(i,3),(cumulant_seeds(i,4)/cumulant_seeds(i,3)^3).^0.5];
        end
            
end
