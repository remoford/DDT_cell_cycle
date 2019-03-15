function [y2] = f_compare(x,m1,s1,m2,s2,bin,h,EType)
    L1=mexWrap('twostage',x,m1,s2,m2,s2,'foo',0.1,'bar');
            
    y2=convolv_2invG_Dirac_option(x,m1,s1,m2,s2,bin,h,EType);
    L2=sum(log(y1));
    
    error=abs(L1-L2);
    fprintf('%f\n',error);
end

