function [pd_max,max_ld]=IMT_analysis_April2017(argstruct)
% This file fits EMG, one-, two-, and three-stage stochasttic models to IMT data.

%argstruct.data is the name of the data, argstruct_model is the model we
%want to fit.
startIMT_analysis=tic;

%This section of code gets the IMT data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(argstruct.data,'erlot')
%for erlotinib
load('experimental_data/erlot_imts_April2017.mat')
data=imt_b;
end

if strcmp(argstruct.data,'at1')
%for AT1 
load('experimental_data/AT1_imts_April2017.mat')
data=imt_b;
end

if strcmp(argstruct.data,'mcf')
%for MCF
load('experimental_data/MCF_imts_April2017.mat')
data=imt_b;
end

if strcmp(argstruct.data,'dmso')
%for DMSO 
load('experimental_data/DMSO_imts_April2017.mat')
data=imt_b;
end

if strcmp(argstruct.data,'chx')
%for CHX
load('experimental_data/CHX_imts_April2017.mat')
data=imt_b;
end

%for FUCCI data
%GetProcessedDataParts
if strcmp(argstruct.data,'fucci')
load('experimental_data/FUCCI_April2017.mat')
data=imt_b;
%data=G2Time_b;
%data=G1Time_b;
end

if strcmp(argstruct.data,'pc9')
%for PC9 cells
load('experimental_data/PC9_April2017.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%choose model to fit
if strcmp(argstruct.model,'two')
[pd_max,max_ld]=twostagefit(data);
end
if strcmp(argstruct.model,'one_lag')
[pd_max,max_ld]=onestagelagfit(data);
end
if strcmp(argstruct.model,'one')
[pd_max,max_ld]=onestagefit(data);
end
if strcmp(argstruct.model,'three')
[pd_max,max_ld]=threestagefit(data);
end
if strcmp(argstruct.model,'emg')
[pd_max,max_ld]=emgfit(data);
end
if strcmp(argstruct.model,'two_lag')
[pd_max,max_ld]=twostagelagfit(data);
end
if strcmp(argstruct.model,'noreset')
[pd_max,max_ld]=twostagefitnoresetfit(data);
end

fprintf("Total runtime:\n")
toc(startIMT_analysis)

end











