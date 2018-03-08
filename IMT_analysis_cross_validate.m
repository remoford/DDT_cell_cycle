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


% GO FAST!
%TolFun = 100;
%TolX = 10;

% GO SLOW!
TolFun = 1;
TolX = 0.001;


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

%fprintf('Selected data:\n');
%data

%select data for the purpose of training and data for the purpose of cross
%validating

dataperm=randperm(length(data));
%size of training data set
trainsize=3*floor(length(data)/4);

%take the first 3/4 of the permuted data for training
datatrain=data(dataperm(1:trainsize));

%keep the remaining data for cross validation
datacross=data(dataperm(trainsize+1:length(dataperm)));

%optimize noreset model
[pd_max_noreset, pd_maxflag_noreset, lcross_noreset, lcrossflag_noreset] = noreset_optimize(datatrain, datacross, TolFun, TolX)

%optimize twostage model
[pd_max, pd_maxflag, lcross, lcrossflag] = twostage_optimize(datatrain, datacross, TolFun, TolX)

%rel is the relative probability of a model compared to another
[AICc rel]=akaikec([lcross lcross_noreset],length(datacross),[4 5])

%make some plots
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

% save the results to a file
save(strcat('crossvalid_results/',dataset,'.mat'));

fprintf("Total runtime:\n")
toc(startIMT_analysis)

end


