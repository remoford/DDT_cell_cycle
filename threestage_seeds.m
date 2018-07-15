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
        vry1 = [.01 .5];
        % proportions of moments represented by the remaining parts
        vry2 = 1-vry1;
        
        c1 = C1*vry;
        c2 = C2*vry;
        c1_comp=C1*vry2;
        c2_comp=C2*vry2;
        m1 = 1./c1;
        m2= 1./c1_comp;
        
%We assume the 2 stage code tries the same number of values as this code
%for each parameter.
        N = length(vry);
        
        P=zeros(N^4,6);
        
        for i=1:N
            for j=1:N
            [P2]=twostage_seeds(c1_comp(i),c2_comp(j));
            for k=1:N^2
            P(N*(i-1)+j+k,:)=[m1(i),s1(j),P2(k,1),P2(k,2),P2(k,3),P2(k,4)];
            end
            end
        end
            
end
