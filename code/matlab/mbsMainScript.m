%% This script calls all lines of code that were involved in the data 
%% analysis of Katyal et al
%%

%% add paths - modify baseFolder

clear

if ispc
    baseFolder = 'C:\Users\skatyal\OneDrive - University College London';
else
    baseFolder = '~/OneDrive - University College London';
end

studyFolder = fullfile(baseFolder, "Projects/Experiments/metaBiasShift/");

addpath(genpath(fullfile(studyFolder, 'analysis')))

cd(studyFolder)

%%
%% Experiment 1 
%%
%% load the raw prolific data into a data structure

mbsDataJsonExp1 = mbsLoadData('exp1_pre', studyFolder);
mbsDataJsonExp1 = mbsLoadData('exp1_full', studyFolder, mbsDataJsonExp1);
mbsDataJsonExp1 = mbsLoadData('exp1_remaining', studyFolder, mbsDataJsonExp1);
mbsDataJsonExp1 = mbsLoadData('exp1_remaining2', studyFolder, mbsDataJsonExp1);

paramsExp1 = mbsSetParams;

%% extract the relevant data from data structure 

exclusionFor = 'mainAnalysis';

mbsDataExp1 = mbsExtractSubj_Exp1(mbsDataJsonExp1, paramsExp1, exclusionFor);

%% export to csv file for running regression models and figures in R

vars2Save = {'conf', 'confZ','corr', 'rt1', 'rt2', ...
    'feedback', 'incdec', 'spe', 'invalid_rt1'};
nLags = 0;

rFilename = 'mbsExp1.csv';

mbsExport2R_Exp1(mbsDataExp1, paramsExp1, vars2Save, rFilename, nLags);

%% fit the different learning models

% models 0-4 : 5 possible non-AD-distortion models (see Supplementary
% Material) - winning model: 4
% models 5-7 : model 4 w feedback, confidence and additive distortions
% models 8-10: model 4 w fb+conf, fb+add, conf+add dist

modelName = 'model_6'; 

nSamples = 10;
zconf = true; % z-score confidence
mhqNum = 2; % 1-phq, 2-gad, 3-spin

mbsDataExp1 = mbsFitSPE(mbsDataExp1, modelName, nSamples, ...
    zconf, mhqNum);

%% save the data structure
save(fullfile(studyFolder, 'data/exp1/processed/mbsDataExp1'), ...
    'mbsDataExp1')

%% compute meta-d' for each subject and task

% add path for HMeta-d toolbox
if ispc
    addpath('C:\Users\skatyal\OneDrive - University College London\Projects\NeuroMatlabToolboxes/HMeta-d/Matlab/')
else
    addpath('/Users/skatyal/OneDrive - University College London/Projects/NeuroMatlabToolboxes/HMeta-d/Matlab/')
end

load(fullfile(studyFolder, 'data/exp1/processed/mbsDataExp1'), 'mbsDataExp1')

mbsDataExp1 = mbsRunMetaDBaseline(mbsDataExp1);

save(fullfile(studyFolder, 'data/exp1/processed/mbsDataExp1'), 'mbsDataExp1')

%%
%%
%% Experiment 2 
%%
%% load data

mbsDataJsonExp2 = mbsLoadData('exp2_pre', studyFolder);
mbsDataJsonExp2 = mbsLoadData('exp2_full', studyFolder, mbsDataJsonExp2);

paramsExp2 = mbsSetParams;
    
%% extract data

exclusionFor = 'mainAnalysis';

[mbsDataExp2, sretChoices] = mbsExtractSubj_Exp2(...
    mbsDataJsonExp2, paramsExp2, exclusionFor, studyFolder);

%% compute factor scores using code from Hopkins et al in R and then read them into the data structure

factorScores = readtable('data/exp2/processed/factor_scores.csv');
mbsDataExp2.questAvg = table2array( factorScores(1:end,2:4));

%% export to R

vars2Save = {'conf', 'confZ', 'corr', 'rt1', 'rt2', ...
    'feedback', 'incdec', 'spe', 'invalid_rt1'};

rFilename = 'mbsExp2.csv';

mbsExport2R_Exp2(mbsDataExp2, paramsExp2, vars2Save, rFilename);

%% fit different learning models

% models 0-4 : 5 possible non-AD-distortion models (see Supplementary
% Material) - winning model: 4
% models 5-7 : model 4 w feedback, confidence and additive distortions
% models 8-10: model 4 w fb+conf, fb+add, conf+add dist
% model 11   : model 4 w conf+add distortion along the 3 transdiag axes

modelName = 'model_6'; 

nSamples = 10;
zconf = true; % if to zscore confidence within participants
mhqNum = 1; % 1-AD, 2-CIT, 3-SW

mbsDataExp2 = mbsFitSPE(mbsDataExp2, modelName, nSamples, ...
    zconf, mhqNum);

%% save data structure with the fits

save(fullfile(studyFolder, 'data/exp2/processed/mbsDataExp2'), 'mbsDataExp2')

%% compute meta-d' for each subject and task

load(fullfile(studyFolder, 'data/exp2/processed/mbsDataExp2'), 'mbsDataExp2')

mbsDataExp2 = mbsRunMetaDBaseline(mbsDataExp2);

save(fullfile(studyFolder, 'data/exp2/processed/mbsDataExp2'), 'mbsDataExp2')
