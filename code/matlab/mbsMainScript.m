%% Script to fit the SPE learning models for Exp 1 and Exp 2
%% Katyal, Huys, Dolan, Fleming - Formation of global confidence and its distortion in anxious-depression
%%

%% add paths - modify baseFolder

clear

addpath(genpath('code'))

%%
%% Experiment 1 
%% load the mat file with data for Exp 1

load('data/exp1/mbsDataExp1.mat')

%% fit the different learning models

% models 0-4 : 5 possible non-AD-distortion models (see Supplementary
% Material) - winning model: 4
% models 5-7 : model 4 w feedback, confidence and additive distortions
% models 8-10: model 4 w fb+conf, fb+add, conf+add dist

modelName = 'model_6'; 

nSamples = 2000;
zconf = true; % z-score confidence
mhqNum = 3; % 1-phq, 2-gad, 3-spin

mbsDataExp1 = mbsFitSPE(mbsDataExp1, modelName, nSamples, ...
    zconf, mhqNum);


%%
%%
%% Experiment 2 
%%
%% load the mat file with data for Exp 1

load('data/exp1/mbsDataExp2.mat')

%% fit different learning models

% models 0-4 : 5 possible non-AD-distortion models (see Supplementary
% Material) - winning model: 4
% models 5-7 : model 4 w feedback, confidence and additive distortions
% models 8-10: model 4 w fb+conf, fb+add, conf+add dist
% model 11   : model 4 w conf+add distortion along the 3 transdiag axes

modelName = 'model_6'; 

nSamples = 2000;
zconf = true; % if to zscore confidence within participants
mhqNum = 1; % 1-AD, 2-CIT, 3-SW

mbsDataExp2 = mbsFitSPE(mbsDataExp2, modelName, nSamples, ...
    zconf, mhqNum);

