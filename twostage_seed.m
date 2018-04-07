function [pd_max, pd_maxflag, lcross, lcrossflag ] = twostage_optimize(datatrain, datacross, TolFun, TolX)
%get sample statistics for fitting initializing the model parameters
num = length(datatrain);
C1 = mean(datatrain);
C2 = var(datatrain);
C3 = sum((datatrain-C1).^3)/(length(datatrain));

numseeds=50

% optimize parameters
pd=zeros(numseeds,4);
ld = NaN*ones(numseeds,1);
flag=zeros(numseeds,1);




for i=1:numseeds
    startOptimization=tic;
    %x0 = P(i,:);
    
    %x0(1) = 0.05225 +(rand(1)-0.5)*2*0.00005;
    %x0(2) = 0.028+(rand(1)-0.5)*2*0.005;
    
    mratio = rand(1);
    sratio = rand(1);
    x0(1) = 1/(C1*mratio);
    x0(3) = 1/(C1*(1-mratio));
    %x0(1) = rand(1);
    %x0(3) = rand(1);
    
    %x0
    %(x0(2)^(2/3))*x0(1)
    x0(2) = rand(1);
    x0(4) = rand(1);
    
    fprintf("optimizing seed %d: m1=%f s1=%f m2=%f s2=%f ", i, x0(1),x0(2),x0(3),x0(4));
    fprintf("ll=%f\n",sum(log(convolv_2invG_adapt_nov(datatrain,x0(1),x0(2),x0(3),x0(4),0.1))));
    
    f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_nov(x,m1,s1,m2,s2,.01);
    myll=@(params)loglikelihood(datatrain, f, 4, params);
    objfun=@(params)penalize(myll, 4, params, [realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax]);
    
    fminsearch_options = optimset('TolFun',TolFun, 'TolX', TolX);
    p=fminsearch(objfun,x0,fminsearch_options);
    
    fprintf("optimized: m1=%f s1=%f m2=%f s2=%f ", p(1),p(2),p(3),p(4));
    
    pd(i,:)=p;
    [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_nov(datatrain,p(1),p(2),p(3),p(4),.01);
    l=sum(log(l));
    
    fprintf("ll=%f\n",l);
    
    ld(i)=l;
    
    m1log(i)=1/p(1);
    s1log(i)=p(2)/(p(1))^(3/2);
    m2log(i)=1/p(3);
    s2log(i)=p(4)/(p(3))^(3/2);
    
    %m1log(i)=1/m1log(i);
    %s1log(i)=(s1log(i)^(2/3))*m1log(i);
    %m2log(i)=1/m2log(i);
    %s2log(i)=(s2log(i)^(2/3))*m2log(i);
    
    lllog(i)=sum(log(convolv_2invG_adapt_nov(datatrain,p(1),p(2),p(3),p(4),0.1)));
    
    fprintf("mean sum delta %f", m1log(i) + m2log(i) - C1);
    fprintf("  sd sum       %f\n", s1log(i) + s2log(i));
    
    if s1log(i) > s2log(i)
        tmp = m1log(i);
        m1log(i) = m2log(i);
        m2log(i) = tmp;
        tmp = s1log(i);
        s1log(i) = s2log(i);
        s2log(i) = tmp;
        tmp = pd(i,1);
        pd(i,1) = pd(i,3);
        pd(i,3) = tmp;
        tmp = pd(i,2);
        pd(i,2) = pd(i,4);
        pd(i,4) = pd(i,2);
    end
    
    toc(startOptimization)
    fprintf("\n");
end

%[foo, I] = sort(lllog);

%best=I(45:end)

%m1log = m1log(best)
%s1log = s1log(best);
%m2log = m2log(best);
%s2log = s2log(best);
%pdn(:,1)=pd(best,1);
%pdn(:,2)=pd(best,2);
%pdn(:,3)=pd(best,3);
%pdn(:,4)=pd(best,4);
%lllog=lllog(best)

figure
hold off
scatter3(m1log,s1log,lllog,'blue');
hold on
scatter3(m2log,s2log,lllog,'red');
%hold on
%scatter3(m1log(best),s1log(best),lllog(best),'green');
%hold on
%scatter3(m2log(best),s2log(best),lllog(best),'green');
title('mean & stddev');

figure
hold off
scatter3(pd(:,1),pd(:,2),lllog,'blue');
hold on
scatter3(pd(:,3),pd(:,4),lllog,'red');
%hold on
%scatter3(pdn(:,1),pdn(:,2),lllog(best),'green');
%hold on
%scatter3(pdn(:,3),pdn(:,4),lllog(best),'green');
title('m & s');

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