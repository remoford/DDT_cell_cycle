% This codes processes IMT data according to cell cycle parts as determined by
% FUCCI and removes cells that are lost and questionable.
% It also selects a cuttoff for the birth time, in order to 
% minimize the bias associated with censored data, due to cells reaching the 
% end of the movie without finishing dividing.

hold off

load('C10_G1_SG2M.mat')

%----------------------Process data----------------------------------------

%remove cells that are not associated with a birthtime
born=find(isnan(Birthtimeh)~=1);

Birthtimeh=Birthtimeh(born);
Lifetimeh=Lifetimeh(born);
G1_time=G1_time(born);
Questionable1=Questionable1(born);
Questionable2=Questionable2(born);
died=died(born);
EOM=EOM(born);
FUCCI_ON_FRAME=FUCCI_ON_FRAME(born);
S_G2_M_time=S_G2_M_time(born);

%remove questionable cells
Q=find(strcmp(Questionable2, {''}));

Birthtimeh=Birthtimeh(Q);
Lifetimeh=Lifetimeh(Q);
G1_time=G1_time(Q);
died=died(Q);
EOM=EOM(Q);
Questionable1=Questionable1(Q);
FUCCI_ON_FRAME=FUCCI_ON_FRAME(Q);
S_G2_M_time=S_G2_M_time(Q);

%remove questionable cells
Q=find(strcmp(Questionable1, {''}));

Birthtimeh=Birthtimeh(Q);
Lifetimeh=Lifetimeh(Q);
G1_time=G1_time(Q);
died=died(Q);
EOM=EOM(Q);
FUCCI_ON_FRAME=FUCCI_ON_FRAME(Q);
S_G2_M_time=S_G2_M_time(Q);

%remove cells that die
dd=find(strcmp(died, 'TRUE')~=1);

Birthtimeh=Birthtimeh(dd);
Lifetimeh=Lifetimeh(dd);
G1_time=G1_time(dd);
died=died(dd);
EOM=EOM(dd);
FUCCI_ON_FRAME=FUCCI_ON_FRAME(dd);
S_G2_M_time=S_G2_M_time(dd);

%segregate data for cells that did and did not reach EOE without FUCCI turning on.
eoe=find(strcmp(EOM, {'TRUE'}));
neoe=find(strcmp(EOM, {'FALSE'}));

%vector of times at which S_G2_M starts for cells that finish S_G2_M
start=Birthtimeh(neoe);
imt=Lifetimeh(neoe);
G2Time=S_G2_M_time(neoe);
G1Time=G1_time(neoe);

%vector of times at which S_G2_M starts for cells that do not finish S_G2_M
start_eoe=Birthtimeh(eoe);
imt_eoe=Lifetimeh(eoe);
G2Time_eoe=S_G2_M_time(eoe);
G1Time_eoe=G1_time(eoe);
%---------------------------------------------------------------------------

%Organize data according to start time

start_time=(min(start)+.1:.1:max(start));

%------------------------------------------------------------------------
% This part of the code finds the last birth time when the proportion of cells 
% born before that time that reach the end of the experiment 
% without dividing is small (cutoff2)

 c1=zeros(length(start_time),1);
 c2=zeros(length(start_time),1);
        for i=1:length(start_time)
            c1(i)=sum(start<=start_time(i));
            c2(i)=sum(start_eoe<=start_time(i));
        end
   %for each start time get proportion of cells born before that time 
   %that reach the end of the experiment without dividing.     
   p_eoe=c2./(c1+c2);
   
  cutoff2=find(p_eoe<=.03, 1, 'last');
 
if isempty(cutoff2)~=1 
    cutoff2=start_time(cutoff2);
    cut2=find(start<=cutoff2);
   
    start_b=start(cut2);
    G2Time_b=G2Time(cut2);
    G1Time_b=G1Time(cut2);
    imt_b=imt(cut2);
    
    else
        G2Time_b=[];
end

%Find the last time that the birthtime and IMT are not
%significantly correlated (cutoff1).  Only consider cells born before
%cutoff1.
if isempty(G2Time_b)~=1
    start_time=(min(start_b)+.1:.1:cutoff2);
    Pvals=zeros(length(start_time),1);
    Corrcoeffs=zeros(length(start_time),1);
    
for i=1:length(start_time)
    ind_i=find(start_b<=start_time(i));
    imt_i=G2Time_b(ind_i);
    start_i=start_b(ind_i);
    if length(imt_i)>2
        
    data_mean=[start_i imt_i];

    [Rm,Pm]=corr(start_i,imt_i,'type','Spearman');
    Pvals(i)=Pm;
    Corrcoeffs(i)=Rm;
    
    RHO = array2table([Rm, Pm], ...
    'VariableNames',{'rho','pval'});
    disp('Correlation between IMT and start time')
    disp(RHO)
   
    end
    
end
cutoff1=find(Pvals>=.01,1,'last')
cutoff1=start_time(cutoff1);
cut1=find(start_b<=cutoff1);
   
    start_b=start_b(cut1);
    G2Time_b=G2Time_b(cut1);
    G1Time_b=G1Time_b(cut1);
    imt_b=imt_b(cut1);

plot(start(start<=cutoff1),G2Time(start<=cutoff1),'bo', 'MarkerSize', 10);
hold on
plot(start(start>cutoff1),G2Time(start>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
plot(start_eoe(start_eoe<=cutoff1),G2Time_eoe(start_eoe<=cutoff1),'go', 'MarkerSize', 10);
plot(start_eoe(start_eoe>cutoff1),G2Time_eoe(start_eoe>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
xlabel('Start Time (h)','FontSize', 20);
ylabel('S-G2-M time (h)','FontSize', 20);
set(gca,'FontSize',20);
saveas(gcf,'G2plot.fig')

hold off

plot(start(start<=cutoff1),G1Time(start<=cutoff1),'bo', 'MarkerSize', 10);
hold on
plot(start(start>cutoff1),G1Time(start>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
plot(start_eoe(start_eoe<=cutoff1),G1Time_eoe(start_eoe<=cutoff1),'go', 'MarkerSize', 10);
plot(start_eoe(start_eoe>cutoff1),G1Time_eoe(start_eoe>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
xlabel('Birth Time (h)','FontSize', 20);
ylabel('G1 time (h)','FontSize', 20);
set(gca,'FontSize',20);
saveas(gcf,'G1plot.fig')


hold off

plot(start(start<=cutoff1),imt(start<=cutoff1),'bo', 'MarkerSize', 10);
hold on
plot(start(start>cutoff1),imt(start>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
plot(start_eoe(start_eoe<=cutoff1),imt_eoe(start_eoe<=cutoff1),'go', 'MarkerSize', 10);
plot(start_eoe(start_eoe>cutoff1),imt_eoe(start_eoe>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
xlabel('Birth Time (h)','FontSize', 20);
ylabel('IMT (h)','FontSize', 20);
set(gca,'FontSize',20);
saveas(gcf,'imtplot.fig')
end