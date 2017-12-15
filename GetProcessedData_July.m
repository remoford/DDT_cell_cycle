% this codes processes IMT data for MCF cells treated with
% CHX and DMSO, and PC9 cells in order to remove cells that are lost and 
% minimize the bias associated with censored data.

hold off

%load('C04.mat')
load('C05.mat')
%load('CHXdata_July.mat')
%load('DMSOdata_July.mat')

%----------------------Process data----------------------------------------

%remove data for cells that leave frame
inframe=find(strcmp(Leftframe, {''}));

Birthtimeh=Birthtimeh(inframe);
Lifetimeh=Lifetimeh(inframe);
EOM=EOM(inframe);
Questionable=Questionable(inframe);
dies=dies(inframe);

%remove cells that are not associated with a birthtime
born=find(isnan(Birthtimeh)~=1);

Birthtimeh=Birthtimeh(born);
Lifetimeh=Lifetimeh(born);
Questionable=Questionable(born);
dies=dies(born);
EOM=EOM(born);

%remove cells that are not associated with a lifetime
divide=find(isnan(Lifetimeh)~=1);

Birthtimeh=Birthtimeh(divide);
Lifetimeh=Lifetimeh(divide);
Questionable=Questionable(divide);
dies=dies(divide);
EOM=EOM(divide);

%remove questionable cells
Q=find(strcmp(Questionable, {''}));

Birthtimeh=Birthtimeh(Q);
Lifetimeh=Lifetimeh(Q);
dies=dies(Q);
EOM=EOM(Q);

%remove cells that die
dd=find(strcmp(dies, 'Y')~=1);

Birthtimeh=Birthtimeh(dd);
Lifetimeh=Lifetimeh(dd);
dies=dies(dd);
EOM=EOM(dd);

%segregate data for cells that did and did not reach EOE without dividing.

neoe = find(strcmp(EOM, 'Y')==0);

eoe=find(strcmp(EOM, 'Y'));

start=Birthtimeh(neoe);
imt=Lifetimeh(neoe);

start_eoe=Birthtimeh(eoe);
imt_eoe=Lifetimeh(eoe);


%---------------------------------------------------------------------------

%compute partial correlation between imt and start time

data=[start imt];
[rho, pval] = partialcorr(data,'type','Spearman');
%[rho, pval] = partialcorr(data,'type','Pearson');

start_sig=0;
if pval(1,2)<=.01
    start_sig=1;
end

rho = array2table(rho, ...
    'VariableNames',{'start','imt'},...
    'RowNames',{,'start','imt'});

disp('Partial Correlation Coefficients')
disp(rho)

pval = array2table(pval, ...
    'VariableNames',{'start','imt'},...
    'RowNames',{'start','imt'});

disp('P Values')
disp(pval)

%-------------------------------------------------------------------------

%---------------------------------------------------------------------------

%Organize data according to start time

start_time=(min(start)+.1:.1:max(start));

%------------------------------------------------------------------------
% This part of the code finds the last time when the proportion of cells 
% that reach the end of the experiment without dividing is small (cutoff2)

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
    imt_b=imt(cut2);
    else
        imt_b=[];
end

%Find the last time that the birthtime and lifetime are not
%significantly correlated (cutoff2).  Only consider cells born before
%cutoff1.
if isempty(imt_b)~=1
    start_time=(min(start_b)+.1:.1:cutoff2);
    Pvals=zeros(length(start_time),1);
    Corrcoeffs=zeros(length(start_time),1);
for i=1:length(start_time)
    ind_i=find(start_b<=start_time(i));
    imt_i=imt_b(ind_i);
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
cutoff1=find(Pvals>=.01,1,'last');
cutoff1=start_time(cutoff1);
cut1=find(start_b<=cutoff1);
   
    start_b=start_b(cut1);
    imt_b=imt_b(cut1);

end

plot(start(start<=cutoff1),imt(start<=cutoff1),'bo', 'MarkerSize', 10);
hold on
plot(start(start>cutoff1),imt(start>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
plot(start_eoe(start_eoe<=cutoff1),imt_eoe(start_eoe<=cutoff1),'go', 'MarkerSize', 10);
plot(start_eoe(start_eoe>cutoff1),imt_eoe(start_eoe>cutoff1),'o', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
xlabel('Birth Time (h)','FontSize', 20);
ylabel('IMT (h)','FontSize', 20);
set(gca,'FontSize',20);
xlim([0 40]);
ylim([0 70]);
saveas(gcf,sprintf('DMSO_scatter_2017.fig'));




