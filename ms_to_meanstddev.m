function [mean1 stddev1 mean2 stddev2] = ms_to_meanstddev(m1, s1, m2, s2)

mean1 = 1/m1;
mean2 = 1/m2;

stddev1 = s1/(m1^(3/2));
stddev2 = s2/(m2^(3/2));


end



