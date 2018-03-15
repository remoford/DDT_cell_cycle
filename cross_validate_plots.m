function [mean1_part1, mean1_part2, sd1_part1, sd1_part2, mean2_part1, mean2_part2, sd2_part1, sd2_part2, rel_avg]=cross_validate_plots(p1,p2,rel)
%p1=cell of parameters for the simple model
%p2=cell of paremeters for the no reset model
%rel=cell of relative probabilities of the models

%meani_partj is a vector giving the mean times spent in part j of the ith
%model

%sdi_partj is a vector giving the standard deviation of times spent in 
%part j of the ith model

% i=1 is the simple model, i=2 is the no reset model

%rel_avg gives the average relative probabilities over all data cuts.

rel_matrix=zeros(10,2);
for i=1:10 
    rel_matrix(i,:)=rel{i};
end
rel_avg=sum(rel_matrix)/10;


p1_matrix=zeros(10,4);
for i=1:10 
    p1_matrix(i,:)=p1{i};
end

p2_matrix=zeros(10,5);
for i=1:10 
    p2_matrix(i,:)=p2{i};
end

mean1_part1=zeros(10,1);
for i=1:10 
    mean1_part1(i,:)=1/p1_matrix(i,1);
end

mean1_part2=zeros(10,1);
for i=1:10 
    mean1_part2(i,:)=1/p1_matrix(i,3);
end

mean2_part1=zeros(10,1);
for i=1:10 
    mean2_part1(i,:)=1/p2_matrix(i,1);
end

mean2_part2=zeros(10,1);
for i=1:10 
    mean2_part2(i,:)=1/p2_matrix(i,3);
end

sd1_part1=zeros(10,1);
for i=1:10 
    sd1_part1(i,:)=p1_matrix(i,2)/(p1_matrix(i,1))^(3/2);
end

sd1_part2=zeros(10,1);
for i=1:10 
    sd1_part2(i,:)=p1_matrix(i,4)/(p1_matrix(i,3))^(3/2);
end

sd2_part1=zeros(10,1);
for i=1:10 
    sd2_part1(i,:)=p2_matrix(i,2)/(p2_matrix(i,1))^(3/2);
end

sd2_part2=zeros(10,1);
for i=1:10 
    sd2_part2(i,:)=p2_matrix(i,4)/(p2_matrix(i,3))^(3/2);
end

end
