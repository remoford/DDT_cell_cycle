function [m1 s1 m2 s2] = meanstddev_to_ms(mean1, stddev1, mean2, stddev2)

m1 = 1/mean1;
m2 = 1/mean2;

s1 = stddev1*(m1^(2/3));
s2 = stddev2*(m2^(2/3));

end