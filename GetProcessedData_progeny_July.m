% this codes asseses the correlation between generation and start time and
% imt and then gets IMT data for cells born before these
% correlations become significant.  It is a condensed version of
% imt_analysis_Jun25 in that it does not analyze the data further by
% generation. 
hold off
%load data
load('AT1.mat')
%load('MCF.mat')
%----------------------Process data----------------------------------------

%remove data for cells that leave frame
inframe=find(LeftFrame==0);

Generation=Generation(inframe);
Parent2=Parent2(inframe);
ID2=ID2(inframe);
CSplit=CSplit(inframe);
CStart=CStart(inframe);
CLifespan=CLifespan(inframe);
Valid=Valid(inframe);

%remove data from F cells
Tcell = find(strcmp(Valid, 'T'));

Generation=Generation(Tcell);
Parent2=Parent2(Tcell);
ID2=ID2(Tcell);
CSplit=CSplit(Tcell);
CStart=CStart(Tcell);
CLifespan=CLifespan(Tcell);
Valid=Valid(Tcell);

genzero = find(strcmp(Parent2, '0'));
Generation(genzero)=600;
Generation(Generation~=600)=NaN;
Generation(genzero)=0;



%insert missing generations
for i=0:10
    parent_ind=find(Generation==i);
    parent_ids=ID2(parent_ind);
    for j=1:length(Parent2)
        flag=find(strcmp(parent_ids,Parent2{j}));
        if isempty(flag)~=1
            Generation(j)=i+1;
        end
    end
end

%remove offspring of invalid cells

legit=find(isnan(Generation)==0);

Generation=Generation(legit);
Parent2=Parent2(legit);
ID2=ID2(legit);
CSplit=CSplit(legit);
CStart=CStart(legit);
CLifespan=CLifespan(legit);

%for MCF divide start by 60
%for AT1 divide start by 10
plot(Generation, CLifespan, 'o');
saveas(gcf,'gene');
plot(CStart./10, CLifespan, 'o');
saveas(gcf,'starte');


%segregate data for cells that did and did not reach EOE without dividing.

eoe=max(CSplit)-5;
eoe_ind=find(CSplit>=eoe);
neoe_ind=find(CSplit<eoe);


Gen=Generation(neoe_ind);
Gen_eoe=Generation(eoe_ind);

Parent=Parent2(neoe_ind);
Parent_eoe=Parent2(eoe_ind);

ID=ID2(neoe_ind);
ID_eoe=ID2(eoe_ind);

Split=CSplit(neoe_ind);
Split_eoe=CSplit(eoe_ind);

Start=CStart(neoe_ind);
Start_eoe=CStart(eoe_ind);

Lifespan=CLifespan(neoe_ind);
Lifespan_eoe=CLifespan(eoe_ind);

%remove cells from generation zero

genzero=find(Gen~=0);
genzero_eoe=find(Gen_eoe~=0);
gen=Gen(genzero);
parent=Parent(genzero);
id=ID(genzero);
split=Split(genzero);
start=Start(genzero);
imt=Lifespan(genzero);

gen_eoe=Gen_eoe(genzero_eoe);
parent_eoe=Parent_eoe(genzero_eoe);
id_eoe=ID_eoe(genzero_eoe);
split_eoe=Split_eoe(genzero_eoe);
start_eoe=Start_eoe(genzero_eoe);
imt_eoe=Lifespan_eoe(genzero_eoe);

%divide by 60 for MCF and 10 for AT1
plot(gen, imt, 'o');
saveas(gcf,'gen');
plot(start./10, imt, 'o');
saveas(gcf,'start');


%---------------------------------------------------------------------------

%compute partial correlation between imt, gen, and start time

variables=[gen start imt];
[rho, pval] = partialcorr(variables,'type','Spearman');
%[rho, pval] = partialcorr(data,'type','Pearson');

gen_sig=0;
start_sig=0;
if pval(3,1)<=.01
    gen_sig=1;
end
if pval(3,2)<=.01
    start_sig=1;
end

rho = array2table(rho, ...
    'VariableNames',{'generation','start','imt'},...
    'RowNames',{'generation','start','imt'});

disp('Partial Correlation Coefficients')
disp(rho)

pval = array2table(pval, ...
    'VariableNames',{'generation','start','imt'},...
    'RowNames',{'generation','start','imt'});

disp('P Values')
disp(pval)

%-------------------------------------------------------------------------

%---------------------------------------------------------------------------

%Organize data according to start time
%step by 6 for MCF and 1 for AT1

start_time=(min(start)+1:1:max(start));
%start_time=(min(start)+6:6:max(start));

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

%We find the last time that the birthtime and lifetime are not
%significantly correlated (cutoff1).  We only consider cells born before
%cutoff2.
if isempty(imt_b)~=1
    % for AT1
    start_time=(min(start_b)+1:1:cutoff2);
    %  for MCF
    %start_time=(min(start_b)+6:6:cutoff2);
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
%for MCF divide start by 60
%for AT1 divide start by 10
plot(start(start<=cutoff1)/10,imt(start<=cutoff1),'bo', 'MarkerSize', 10);
hold on
plot(start(start>cutoff1)/10,imt(start>cutoff1),'x', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
plot(start_eoe(start_eoe<=cutoff1)/10,imt_eoe(start_eoe<=cutoff1),'g+', 'MarkerSize', 10);
plot(start_eoe(start_eoe>cutoff1)/10,imt_eoe(start_eoe>cutoff1),'x', 'MarkerSize', 10, 'Color',[.5,.5,.5]);
xlabel('Birth Time (h)','FontSize', 20);
ylabel('IMT (h)','FontSize', 20);
set(gca,'FontSize',20);
