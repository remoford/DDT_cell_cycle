function cross_validate_plots(p1,p2,rel,data)
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

avgRelSimple = 0;
for i=1:10
   avgRelSimple = avgRelSimple + rel{i}(1);
end
avgRelSimple = avgRelSimple/10;

avgRelMeyer = 0;
for i=1:10
   avgRelMeyer = avgRelMeyer + rel{i}(2);
end
avgRelMeyer = avgRelMeyer/10;

maxMean = max(max(max(mean1_part1),max(mean2_part1)),max(max(mean1_part2),max(mean2_part2)));
maxMeanPadded = maxMean + 0.1*maxMean;
maxStdDev = max(max(max(sd1_part1),max(sd2_part1)),max(max(sd1_part2),max(sd2_part2)));
maxStdDevPadded = maxStdDev + 0.1*maxStdDev;
% close all existing figures
delete(findall(0,'Type','figure'))


figure('Name','Data and Fits');
hold on;
[counts,centers] = hist(data, 20);
tot=sum(counts);
wdth=centers(2)-centers(1);
hght=counts/(tot*wdth);
bar(centers,hght, 'FaceColor', 'white', 'LineWidth', 1.5)
hold on
tt=min(data):.01:max(data);
for kk=1:10
    plot(tt,convolv_2invG_adapt_nov(tt,p1{kk}(1),p1{kk}(2),p1{kk}(3),p1{kk}(4),.01),'b');
    plot(tt,convolv_2invG_noreset(tt,p2{kk}(1),p2{kk}(2),p2{kk}(3),p2{kk}(4),p2{kk}(5),.01),'r');
    xlabel('Time');
    ylabel('Probability');
    title('Data and Fits');
    legend('Data','DDT Model','Meyer Model');
end

colormap jet;
pallete = colormap;
indexSpace = round(linspace(1,64,10));

figure('Name','Parameters Simple');
hold on;
for kk=1:10
    scatter(mean1_part1(kk,:),sd1_part1(kk,:),(rel{kk}(1)+0.1)*200,pallete(indexSpace(kk),:),'LineWidth',1.5); 
    scatter(mean1_part2(kk,:),sd1_part2(kk,:),(rel{kk}(1)+0.1)*200,pallete(indexSpace(kk),:),'LineWidth',1.5);
    xlabel({'Mean';'Note: Color represents the cross validation cut'});
    ylabel('Standard Deviation');
    title({'Drift Diffusion Threshold Model Parameters';sprintf("Average Strength = %f", avgRelSimple)});
    axis([0 maxMeanPadded 0 maxStdDevPadded]);
end

figure('Name','Parameters Meyer');
hold on;
for kk=1:10
    yyaxis left;
    scatter(mean2_part1(kk,:),sd2_part1(kk,:),(rel{kk}(2)+0.1)*200,pallete(indexSpace(kk),:),'LineWidth',1.5); 
    scatter(mean2_part2(kk,:),sd2_part2(kk,:),(rel{kk}(2)+0.1)*200,pallete(indexSpace(kk),:),'LineWidth',1.5);
    ylabel('Standard Deviation');
    title({'Memory (Meyer) Model Parameters';sprintf("Average Strength = %f", avgRelMeyer)});
    axis([0 maxMeanPadded 0 maxStdDevPadded]);
    xlabel({'Mean';'Note: Size represents the relative AIC strength of the model'});
    yyaxis right;
    ylim([0 1]);
    scatter(maxMeanPadded,p2{kk}(5),(rel{kk}(2)+0.1)*200,pallete(indexSpace(kk),:),'LineWidth',1.5);
    ylabel('Probability of Skipping the Restriction Point');
end


end
