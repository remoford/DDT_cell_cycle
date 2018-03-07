function []=IMT_analysis_cross_validate(dataset)

% This file fits two-stage with and without perfect reset models to IMT 
%data and performs cross validation. 75% of the data is used for training,
%the remaining 25% goes to cross validation.

%For each model, as many as two "best" parameter sets are evaluated using
%cross validation, the best choice wherein the distribution of time spent 
%in one part of the cell cycle is a Dirac delta distribution 
%(flagged model) and the best choice where the distribution of time in both 
%parts is inverse Guassias distributed (no flag model).  

%lcross=the log likelihood of the data given the best model with perfect 
%reset and no flag.
%lcrossflag=the log likelihood of the data given the best model with perfect 
%reset and a flag.
%lcross_noreset=the log likelihood of the data given the best model with  
%no reset and no flag.
%lcrossflag_noreset=the log likelihood of the data given the best model with  
%no reset and a flag.

startIMT_analysis=tic;

%This section of code gets the IMT data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch dataset
    case 'erlotinib'
        %for erlotinib
        load('experimental_data/erlot_imts_April2017.mat')
        data=imt_b;
        fprintf('Selected erlotinib dataset\n');

    case 'AT1'
        %for AT1 
        load('experimental_data/AT1_imts_April2017.mat')
        data=imt_b;

    case 'MCF'
        %for MCF
        load('experimental_data/MCF_imts_April2017.mat')
        data=imt_b;

    case 'DMSO'
        %for DMSO 
        load('experimental_data/DMSO_imts_April2017.mat')
        data=imt_b;

    case 'synthetic'
        data=[];
        for i=1:50
            data = [data, random('InverseGaussian',0.3,0.8)];
        end
        data = data';

    case 'CHX'
        %for CHX
        load('experimental_data/CHX_imts_April2017.mat')
        data=imt_b;

    case 'FUCCI'
        %for FUCCI data
        %GetProcessedDataParts
        load('experimental_data/FUCCI_April2017.mat')
        data=imt_b;
        %data=G2Time_b;
        %data=G1Time_b;

    case 'PC9'
        %for PC9 cells
        load('experimental_data/PC9_April2017.mat')
end

fprintf('Selected data:\n');
data

%select data for the purpose of training and data for the purpose of cross
%validating

dataperm=randperm(length(data));
%size of training data set
trainsize=3*floor(length(data)/4);
%take the first 3/4 of the permuted data for training
datatrain=data(dataperm(1:trainsize));
%keep the remaining data for cross validation
datacross=data(dataperm(trainsize+1:length(dataperm)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get sample statistics for fitting initializing the model parameters
num = length(datatrain);
C1 = mean(datatrain);
C2 = var(datatrain);
C3 = sum((datatrain-C1).^3)/(length(datatrain));



[pd_max_noreset, pd_maxflag_noreset, lcross_noreset, lcrossflag_noreset] = noreset_optimize(datatrain, datacross)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fit two-stage model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % BEGIN FUNCTION FIT_TWOSTAGE
    
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
        options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');
        flag=zeros(length(id),1);
        for i=1:length(id)  
            x0 = P(i,:);
            f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_nov(x,m1,s1,m2,s2,.01);
            %f=@(x,m1,s1,m2,s2)convolv_2invG_adapt2(x,m1,s1,m2,s2,.01,4);
            [p,conf1]=mle(datatrain,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options)
            pd(i,:)=p;
            [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_nov(datatrain,p(1),p(2),p(3),p(4),.01);
            l=sum(log(l));
            ld(i)=l    

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
        
        if isempty(indflag)==0
         % common to each fit, consider factoring out
        [max_ldflag,row_ldflag]=max(ld_trueflag)
        %best flagged model
        pd_maxflag = pdflag(row_ldflag,:)
         %perform cross validation
        [lcrossflag]=convolv_2invG_adapt_nov(datacross,pd_maxflag(1),pd_maxflag(2),pd_maxflag(3),pd_maxflag(4),.01);
        lcrossflag=sum(log(lcrossflag));
        end
        
       
        
     
    % END FUNCTION FIT_TWOSTAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rel is the relative probability of a model compared to another
[AICc rel]=akaikec([lcross lcross_noreset],length(datacross),[4 5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[counts,centers] = hist(datacross,20);
tot=sum(counts);
wdth=centers(2)-centers(1);
hght=counts/(tot*wdth);
bar(centers,hght)
hold on
tt=min(data):.01:max(data);
plot(tt,convolv_2invG_adapt_nov(tt,pd_max(1),pd_max(2),pd_max(3),pd_max(4),.01),'b');
plot(tt,convolv_2invG_noreset(tt,pd_max_noreset(1),pd_max_noreset(2),pd_max_noreset(3),pd_max_noreset(4),pd_max_noreset(5),.01),'r');

title('IMT with simple and no reset model')
xlabel('Intermitotic Times (IMT)')
ylabel('Probability density function (pdf)')
legend('cells', 'simple model','no reset model')

fprintf("Total runtime:\n")
toc(startIMT_analysis)

save(strcat('crossvalid_results/',dataset,'.mat'));

end


