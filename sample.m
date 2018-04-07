%for CHX
load('experimental_data/CHX_imts_April2017.mat')
data=imt_b;
for i=1:1000
    m1log(i)=rand(1);
    s1log(i)=rand(1);
    m2log(i)=rand(1);
    s2log(i)=rand(1);
    [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_nov(data,m1log(i),s1log(i),m2log(i),s2log(i),0.1);
    l=sum(log(l));
    lllog(i)=l
end

hold off
scatter3(m1log,s1log,lllog,'blue');
hold on
scatter3(m2log,s2log,lllog,'red');
hold off