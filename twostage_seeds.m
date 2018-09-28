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
        m1 = 1./c1;
        m2= 1./c1_comp;
        s1 = (c2./c1.^3).^0.5;
        s2 = (c2_comp./c1_comp.^3).^0.5;
        N = length(vry1);
        
        P=zeros(N^2,4);
        
        for i=1:N
            for j=1:N
            
            P(N*(i-1)+j,:)=[m1(i),s1(j),m2(i),s2(j)];
            end
        end
            
end
