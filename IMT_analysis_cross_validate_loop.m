function []=IMT_analysis_cross_validate_loop(dataset)

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

%let's edit this some more!

rng(1);

% GO FAST!
%TolFun = 100;
%TolX = 10;

% GO SLOW!
TolFun = 1;
TolX = 0.001;

startIMT_analysis=tic;

varreset_enable=0;
noreset_enable=1;
twostage_enable=1;
analysis_enable=1;

options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs');


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


pd_max_noreset = {};
pd_maxflag_noreset = {};
lcross_noreset = {};
lcrossflag_noreset = {};

pd_max = {};
pd_maxflag = {};
lcross = {};
lcrossflag = {};

AICc = {};
rel = {};
for kk=1:10
    
%select data for the purpose of training and data for the purpose of cross
%validating

dataperm=randperm(length(data));
%size of training data set
trainsize=3*floor(length(data)/4);
%take the first 3/4 of the permuted data for training
datatrain=data(dataperm(1:trainsize));
%keep the remaining data for cross validation
datacross=data(dataperm(trainsize+1:length(dataperm)));
         
if varreset_enable
    varreset_optimize();
end

if noreset_enable
    [pd_max_noreset{kk}, pd_maxflag_noreset{kk}, lcross_noreset{kk}, lcrossflag_noreset{kk}] = noreset_optimize(datatrain, datacross, TolFun, TolX);

end

if twostage_enable
    [pd_max{kk}, pd_maxflag{kk}, lcross{kk}, lcrossflag{kk}] = twostage_optimize(datatrain, datacross, TolFun, TolX);
end

if analysis_enable

    %rel is the relative probability of a model compared to another
    [AICc{kk} rel{kk}]=akaikec([lcross{kk} lcross_noreset{kk}],length(datacross),[4 5])

    [counts,centers] = hist(datacross,20);
    tot=sum(counts);
    wdth=centers(2)-centers(1);
    hght=counts/(tot*wdth);
    bar(centers,hght)
    hold on
    tt=min(data):.01:max(data);
    plot(tt,convolv_2invG_adapt_nov(tt,pd_max{kk}(1),pd_max{kk}(2),pd_max{kk}(3),pd_max{kk}(4),.01),'b');
    plot(tt,convolv_2invG_noreset(tt,pd_max_noreset{kk}(1),pd_max_noreset{kk}(2),pd_max_noreset{kk}(3),pd_max_noreset{kk}(4),pd_max_noreset{kk}(5),.01),'r');
    %plot(tt,convolv_2invG_var_reset(tt,pd_max_var_reset{kk}(1),pd_max_var_reset{kk}(2),pd_max_var_reset{kk}(3),pd_max_var_reset{kk}(4),pd_max_var_reset{kk}(5),.01),'g');

    title('IMT with simple and no reset model')
    xlabel('Intermitotic Times (IMT)')
    ylabel('Probability density function (pdf)')
    legend('cells', 'simple model','no reset model')


    %savefig(gcf,sprintf('IMT_hist%d',kk))

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fprintf("Total runtime:\n")
toc(startIMT_analysis)

save(strcat('crossvalid_results/',dataset,'_loop.mat'));

end


