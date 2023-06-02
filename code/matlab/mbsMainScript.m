%% Script to fit the SPE learning models for Exp 1 and Exp 2
%% Katyal, Huys, Dolan, Fleming - How underconfidence is maintained in anxiety and depression
%%
% BEFORE PROCEEDING: Set current directory to the base directory of
% downloaded repo from github.com/sucharitk/confidence-distortion-AD

%% 1: add paths

clear
close all
addpath(genpath('code'))

%%
%% Experiment 1 
%% 2: load the mat file with data for Exp 1

load('data/exp1/mbsDataExp1.mat')

%% 3: fit the different learning models for Exp 1

% models 0-4 : 5 possible non-AD-distortion models (see Supplementary
% Material) - winning model: 4
% models 5-7 : model 4 w feedback, confidence, and additive distortions
% models 8-10: model 4 w fb+conf, fb+add, conf+add dist

% for running each model from 0-10 modify this string (e.g., model_0,
% model_1... model_10) and execute this cell each time. the results of each
% model will be saved in mbsDataExp1 data structure
modelName = 'model_6'; 

% NOTE: we ran the model with 2000x3 mcmc samples. With this number the
% model will take a long time to run. To speed up, reduce this to a smaller
% value like nSamples = 30 or 50. In general, the means of the estimated
% posteriors are not too affected by this value
nSamples = 2000; 
mhqNum = 1; % 1-phq, 2-gad, 3-spin

mbsDataExp1 = mbsFitSPE(mbsDataExp1, modelName, nSamples, mhqNum);

%% 4: after fitting models 0-10, save the mat file
% the same file will be loaded by mbsAnalyses_model.R to plot the figures

save(fullfile(studyFolder, 'data/exp1/mbsDataExp1'), 'mbsDataExp1')


%%
%%
%% Experiment 2 
%%
%% 5: load the mat file with data for Exp 2

load('data/exp2/mbsDataExp2.mat')

%% 6: fit different learning models

% models 0-4 : 5 possible non-AD-distortion models (see Supplementary
% Material) - winning model: 4
% models 5-7 : model 4 w feedback, confidence and additive distortions
% models 8-10: model 4 w fb+conf, fb+add, conf+add dist
% model 11   : model 4 w conf+add distortion along the 3 transdiag axes

% for running each model from 0-10 modify this string (e.g., model_0,
% model_1... model_10) and execute this cell each time. the results of each
% model will be saved in mbsDataExp1 data structure
modelName = 'model_10'; 

nSamples = 10; 
mhqNum = 1; % 1-AD, 2-CIT, 3-SW

mbsDataExp2 = mbsFitSPE(mbsDataExp2, modelName, nSamples, mhqNum);

%% 7: after fitting models 0-11, save the mat file
% the same file will be loaded by mbsAnalyses_model.R to plot the figures

save(fullfile(studyFolder, 'data/exp2/mbsDataExp2'), 'mbsDataExp2')
