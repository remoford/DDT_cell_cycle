function P=onestagepdf_binned(imt,mu,s,h)
                x=0:h:max(imt);
                C=onestagepdf2(x,mu,s);
                n=length(imt);
                I=zeros(size(C));
                P=zeros(size(imt));

                goback=.1/h;
                
                for i=1:n
                %find element of x that is closest to t(i)
                    [~,I(i)]=min((imt(i)-x).^2);
                    %If t(i)<0 the probability is set to zero, otherwise the
                    %probability is approxiated as a value from the vector x.
                    I_vector=(I(i)-goback+1):1:I(i);
                        
                    if 0 >= min(I_vector)
           
                        I_vector=I_vector(I_vector>0);
            
                    end
                        
                        P(i)=sum(C(I_vector))*h;

                end
                %toc
                %Y=max(realmin,P);
             
end