function [pd_max, pd_maxflag, lcross, lcrossflag] = twostage_optimize(datatrain, datacross, TolFun, TolX)
%get sample statistics for fitting initializing the model parameters
num = length(datatrain);
C1 = mean(datatrain);
C2 = var(datatrain);
C3 = sum((datatrain-C1).^3)/(length(datatrain));

% prepare statistical parameters
vry = [.25 .5 .75]';
% vrys=[.01 1 10]';
% [x0, l_moments, min_res] = moments_method_2stage(data);
% m1 = x0(1)*vry; 
% s1 = x0(2)*vrys;
% m2 = x0(3)*vry;
% s2 = x0(4)*vrys;
c1 = C1*vry;
c2 = C2*vry;
m = 1./c1;
s = (c2./c1.^3).^0.5;
N = length(vry);

% prepare parameter seeds

%get all pairs of the form [m(i),s(j)]
%these pairs represent all possible unique 
%parameter choices for each part of the cell
%cycle.  
pcomb = allcomb(m,s);

%place paramter pairs into a cell.  The parameters choices for each part
%are now indexed
pcell = cell(length(pcomb),1);
for i = 1:length(pcomb)
    pcell{i} = pcomb(i,:);
end

%get all pairs of indices for the parameter 
%choices for each part of the cycle to get all 
%parameter choices for the entire cycle
id = allcomb(1:length(pcomb),1:length(pcomb));

%sort the pairs in ascending order.  
%This equates choices of the form [i,j] and [j,i].
id = sort(id,2);

%remove repeats
id = unique(id,'rows');

%create a matrix of unique parameter choices for the cell cycle
P = zeros(length(id),4);
for ii = 1:length(id)
    P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}];
end

% optimize parameters
pd=zeros(length(P),4);
ld = NaN*ones(length(P),1);
flag=zeros(length(id),1);

for i=1:length(id)  
    startOptimization=tic;
    x0 = P(i,:);

    fprintf("optimizing seed %d: m1=%f s1=%f m2=%f s2=%f\n", i, x0(1),x0(2),x0(3),x0(4));

    f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_window(x,m1,s1,m2,s2,.01);
    %f=@(x,m1,s1,m2,s2)convolv_2invG_adapt2(x,m1,s1,m2,s2,.01,4);
    %[p,conf1]=mle(datatrain,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options)

%%%%%These are the options%%%%%%%%%%%%%%%%%%%
%set max iter to 1
    fminsearch_options = optimset('TolFun',TolFun, 'TolX', TolX);
    myll=@(params)loglikelihood(datatrain, f, 4, params);
    objfun=@(params)penalize(myll, 4, params, [realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax]);
    
    %make a for loop running the optimizer 10000 times so we can see the traejctory.
    for i=1:10000
        p=fminsearch(objfun,x0,fminsearch_options);
        
    %add line to plot p, set hold on so we can see successive p values
    plot(p)
    hold on
    
    %end for loop
    end
    
    %save plot with index indicating value of i (which indentifies the
    %initial data.)

    fprintf("optimized: m1=%f s1=%f m2=%f s2=%f\n", p(1),p(2),p(3),p(4));

    pd(i,:)=p;
    [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_window(datatrain,p(1),p(2),p(3),p(4),.01);
    l=sum(log(l));

    fprintf("log-liklihood=%f\n",l);

    ld(i)=l;

    toc(startOptimization)
    fprintf("\n");
end

% we previously optimized with a larger step size, recalculate with
% a smaller stepsize after the fact
%ld_true=zeros(length(ld),1);
%for i=1:length(ld)
    %[l,hp_true(i),flag_true(i),E_true(i)]=convolv_2invG_adapt_nov(datatrain,pd(i,1),pd(i,2),pd(i,3),pd(i,4),.001);
    %ld_true(i)=sum(log(l));
%end
%find the best flagged model and the best model that is not flagged

%skip recalculation for now
ld_true=ld;
flag_true=flag;

indflag=find(flag_true==1);
indnoflag=find(flag_true==0);

ld_trueflag=ld_true(indflag);
ld_true=ld_true(indnoflag);

pdflag=pd(indflag,:);
pd=pd(indnoflag,:);

% common to each fit, consider factoring out
[max_ld,row_ld]=max(ld_true)

%best nonflagged model
pd_max = pd(row_ld,:)

%perform cross validation
[lcross]=convolv_2invG_adapt_nov(datacross,pd_max(1),pd_max(2),pd_max(3),pd_max(4),.01);
lcross=sum(log(lcross));

pd_maxflag=[];
lcrossflag=[];
if isempty(indflag)==0
    % common to each fit, consider factoring out
    [max_ldflag,row_ldflag]=max(ld_trueflag)
    
    %best flagged model
    pd_maxflag = pdflag(row_ldflag,:)
    
    %perform cross validation
    [lcrossflag]=convolv_2invG_adapt_nov(datacross,pd_maxflag(1),pd_maxflag(2),pd_maxflag(3),pd_maxflag(4),.01);
    lcrossflag=sum(log(lcrossflag));
end

end