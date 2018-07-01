function P = Conv2(z,y,h,n,t,x)

% BEGIN FUNCTION DOTHECONVOLUTION_OUTER
        % Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x
        % Outputs: P
            % BEGIN FUNCTION DOTHECONVOLUTION_INNER
            % Input parameters: z, y, h, n, t, i, I, x
            % Outputs: logP0
                % find the discrete convolution of the vectors y and z
                % the (i-1)th element of v approximates the convolution of the pdfs 
                % over [.001, x(i)] as a left-hand Riemann sum.
                C=conv(z,y)*h;
                N=length(y);
                % only the first N elements of the convolution are valid
                C=C(1:N);
                I=zeros(n,1);
                P=zeros(n,1);
                goback=.1/h;
                for i=1:n
                %find element of x that is closest to t(i)
                    [~,I(i)]=min((t(i)-x).^2);
                    %If t(i)<0 the probability is set to zero, otherwise the
                    %probability is approxiated as a value from the vector x.
                    if t(i)>0 && I(i)>1
                        I_vector=(I(i)-goback+1):1:I(i);
                        
                        if 0 >= min(I_vector)
                            I_vector=I_vector(I_vector>0);
                        end
                        
                        P(i)=sum(C(I_vector))*h;

                    else
                        P(i)=realmin;
                    end
                end
                %toc
                P=max(realmin,P);
            % END FUNCTION DOTHECONVOLUTION_INNER