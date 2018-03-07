function [pd_max_noreset, pd_maxflag_noreset, lcross_noreset, lcrossflag_noreset] = noreset_optimize(datatrain, datacross)


%get sample statistics for fitting initializing the model parameters
num = length(datatrain);
C1 = mean(datatrain);
C2 = var(datatrain);
C3 = sum((datatrain-C1).^3)/(length(datatrain));

%Fit two-stage model with imperfect reset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BEGIN FUNCTION FIT_TWOSTAGE_NORESET

% prepare statistical parameters
vry = [.1 .5 .9]';
r = [.01 .5 .99]';
%c1 and c2 give the possible initial guess means and variances
c1 = C1*vry;
c2 = C2*vry;
%m and s give the possible initial guess parameter values
m = 1./c1;
s = (c2./c1.^3).^0.5;
N = length(vry);

% prepare parameter seeds

%get all pairs of the form [m(i),s(j)]
%these pairs represent all possible unique 
%parameter choices for each part of the cell
%cycle.  
pcomb = allcomb(m,s);
%place paramter triples into a cell.  The parameter choices for each part
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
P = zeros(length(id)*length(r),5);


% STUB CODE STUB CODE
% this needs to be fixed so that there are a variety of starting
% points for the r parameter. As is, this hard codes a starting
% point of 0.5 which is not what we really want!
for ii = 1:length(id)
    for jj=1:length(r)
    P(length(r)*(ii-1)+jj,:) = [pcell{id(ii,1)},pcell{id(ii,2)},r(jj)];
    end
end

% optimal parameters for each initial guess
pd_noreset=zeros(length(P),5);
%likelihoods for each initial guess
ld_noreset = NaN*ones(length(P),1);
options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
%flag is one if a part is approximated as a Dirac delta, so there
%is no varaibility in time spent in that part.
flag_noreset=zeros(length(id),1);
for i=1:(length(id)*length(r))
    startOptimization=tic;
    %sets initial guess
    x0 = P(i,:);
    fprintf("optimizing seed %d: m1=%f s1=%f m2=%f s2=%f r=%f\n", i, x0(1),x0(2),x0(3),x0(4),x0(5));
    f=@(x,m1,s1,m2,s2,r)convolv_2invG_noreset(x,m1,s1,m2,s2,r,.01);
    %[p,conf]=mle(datatrain,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf 1],'lowerbound',[0 0 0 0 0],'options',options);

    fminsearch_options = optimset('TolFun',10, 'TolX', 1);
    myll=@(params)loglikelihood(datatrain, f, 5, params);
    objfun=@(params)penalize(myll, 5, params, [realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax;0.001  0.999]);
    p=fminsearch(objfun,x0,fminsearch_options);


    fprintf("optimized: m1=%f s1=%f m2=%f s2=%f r=%f\n", p(1),p(2),p(3),p(4),p(5));
    %save parameters
    pd_noreset(i,:)=p;
    %gets the likelihhod of the parameters
    [l,hp(i),flag_noreset(i),E(i)]=convolv_2invG_noreset(datatrain,p(1),p(2),p(3),p(4),p(5),.01);
    l=sum(log(l));
    fprintf("log-liklihood=%f\n",l);
    if flag_noreset(i) == 1
        fprintf("used dirac delta approximation in final result\n");
    end
    ld_noreset(i)=l;
    toc(startOptimization)
    fprintf("\n");
end
        
        % we previously optimized with a larger step size, recalculate with
        % a smaller stepsize after the fact
%         ld_true_noreset=zeros(length(ld),1);
%         for i=1:length(ld)
%             [l,hp_true(i),flag_true_noreset(i),E_true(i)]=convolv_2invG_noreset(datatrain,pd_noreset(i,1),pd_noreset(i,2),pd_noreset(i,3),pd_noreset(i,4),pd_noreset(i,5),.001);
%             ld_true_noreset(i)=sum(log(l));
%         end
        
% skip recalculation for now

ld_true_noreset=ld_noreset;
flag_true_noreset=flag_noreset;

%find the best flagged model and the best model that is not flagged

%indices of flagged models
indflag_noreset=find(flag_true_noreset==1);
%ibndices of no flag models
indnoflag_noreset=find(flag_true_noreset==0);


ld_trueflag_noreset=ld_true_noreset(indflag_noreset);
ld_true_noreset=ld_true_noreset(indnoflag_noreset);

pdflag_noreset=pd_noreset(indflag_noreset,:);
pd_noreset=pd_noreset(indnoflag_noreset,:);

% common to each fit, consider factoring out

[max_ld_noreset,row_ld_noreset]=max(ld_true_noreset)
%best nonflagged model
pd_max_noreset = pd_noreset(row_ld_noreset,:)
%cross validate
[lcross_noreset]=convolv_2invG_noreset(datacross,pd_max_noreset(1),pd_max_noreset(2),pd_max_noreset(3),pd_max_noreset(4),pd_max_noreset(5),.01);
lcross_noreset=sum(log(lcross_noreset));

pd_maxflag_noreset=[];
lcrossflag_noreset=[];
if isempty(indflag_noreset)==0
    [max_ldflag_noreset,row_ldflag_noreset]=max(ld_trueflag_noreset)
    %best flagged model
    pd_maxflag_noreset = pdflag_noreset(row_ldflag_noreset,:)
    %perform cross validation
    [lcrossflag_noreset]=convolv_2invG_noreset(datacross,pd_maxflag_noreset(1),pd_maxflag_noreset(2),pd_maxflag_noreset(3),pd_maxflag_noreset(4),pd_maxflag_noreset(5),.01);
    lcrossflag_noreset=sum(log(lcrossflag_noreset));
end
        
        
        
    % END FUNCTION FIT_TWOSTAGE_NORESET

end