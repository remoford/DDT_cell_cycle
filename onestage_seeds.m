function [P]=onestage_seeds(c1,c2)

%This code finds the parameter seed for the onestage model that matches the
%first two moments of the data.

        m = 1/c1;
        
        s = (c2/c1^3)^0.5;
        
        P=[m s];
            
end
