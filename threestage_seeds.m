function [P]=threestage_seeds(C1,C2)

%C1 is the mean of the data
%C2 is the variance of the data

%This code finds parameter seeds for a threestage model that match the mean 
% and variance of the data.  
%These seeds correspond to parameter combinations wherein a single part 
%of the model accounts for 10% or 50% of the variance and 10% or 50% of the
%mean.  The remaining parts of the model account for 90% or 50% of these
%moments, respectively.  

    
        % proportions of moments represented by the first part
        vry1 = [.1 .2 .3 .4 .5];
        vry2 = [.1 .5];
        % proportions of moments represented by the remaining parts
        portion = abs(1-vry1);
        
        c1 = C1*vry1;
        c2 = C2*vry1;
        
        c1_comp=C1*portion;
        c2_comp=C2*portion;
        
        m1 = 1./c1;
        s1 = (c2./c1.^3).^0.5;
        
        
%We assume the 2 stage code tries the same number of values as this code
%for each parameter.
        N1 = length(vry1);
        N2 = length(vry2);
        
        P=zeros((N1^2)*(N2^2),6);
        
        for i=1:N1
            for j=1:N1
            [P2]=twostage_seeds(c1_comp(i),c2_comp(j),vry2);
            for k=1:N2^2
            P((N^3)*(i-1)+(N^2)*(j-1)+k,:)=[m1(i),s1(j),P2(k,1),P2(k,2),P2(k,3),P2(k,4)];
            end
            end
        end
            
end
