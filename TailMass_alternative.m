function check2 = TailMass_alternative(m,s,eps,T2)

% BEGIN FUNCTION TailMass
        % Input parameters: m, s, eps, T2, nu, sd,
        % Outputs: check2
            % to estimate the error (in calculating probabilities the of the
            % data), that results from approximating the first pdf as a
            % point-mass distribtion, find the maximum of the absolute value
            % of the derivative of the second pdf.
            gp=gp_max(m(2),s(2));

            % determine the radius, r, of a small interval over which the
            % 1. second pdf, g, is approximately constant, i.e. changes by less than eps/3
            % over any interval with that radius
            % and 2. g(t) is small for t<r.  ?????
            r=min(eps/(3*gp),T2);
            checkval=onestagepdf2(r,m(2),s(2));
            while checkval>=eps/2
                r=r/2;
                checkval=onestagepdf2(r,m(2),s(2));
            end

            % This is the point at
            % which its mass is concentrated.
            T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1)));

           F1=onestagepdf2(T1-r,m(1),s(1));
           F2=onestagepdf2(T1+r,m(1),s(1));
           
           check2=F1+F2;
            %     end
        % END FUNCTION TailMass