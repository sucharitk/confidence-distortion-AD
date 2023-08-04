function fit = fit_globalSPE(spe,corr,conf,feedback,task, ...
    ntrials,modelNum, nsamples, covariate, tmpfolder, runParallel, doPlot)
%
% model to fit SPE + confidence data
%
%
%
% Variables passed to run the model: 
% 
% spe: end-of-run self-performance estimate [NumSubjects x NumBlocks]
%
% corr: task accuracy [NumSubjects x NumBlocks X NumTrials]
%
% conf: confidence reports [NumSubjects x NumBlocks X NumTrials]
%
% feedback: feedback (+1: positive, 0: none, -1: negative). if no feedback
% in your data, set all to 0
% [NumSubjects x NumBlocks X NumTrials]
%
% task: if using 2 tasks, specify task used on the respective block 
% (accepted values: 1, 2). for 1 task, specify all 1s [NumSubjects x NumBlocks].
%
% ntrials: number of trials on each block (in case some blocks have
% different number of trials than others) [1 x NumBlocks]
%
% modelNum: which model to run (current accepted values: 1-11) specified as
% a string 
%
% nsamples: number of mcmc samples
%
% covariate: individual difference covariate [1 x NumSubjects]
%
% tmpfolder: name of folder to store temporary jags files while running
% model
%
% runParallel: 1-run model using parallel processing, 0-not in parallel
%
% doPlot: plot the posterior distributions (0, 1)
%
%

sz = size(corr);
nSubj = sz(1);
nBlocks = sz(2);

cwd = pwd;

findpath = which('fit_mbs_group_0.txt');
if isempty(findpath)
    error('Please add model directory to the path')
else
    hmmPath = fileparts(findpath);
    cd(hmmPath)
end

if ~exist('tmpfolder', 'var')
    rn = num2str(rand);
    tmpfolder = ['tmpjags' rn(3:end)];
end

if ~exist('doPlot', 'var')
    doPlot = true;
end

newdircreated = false;
% if ~exist(tmpfolder,'dir')
mkdir(tmpfolder);
newdircreated = true;
% end

if ~exist('nsamples','var')
    nsamples = 250;
end

%% Sampling
if ~exist('mcmc_params','var')
    % MCMC Parameters
    mcmc_params.nchains = 3; % How Many Chains?
    mcmc_params.nburnin = nsamples/2; % How Many Burn-in Samples?
    mcmc_params.nsamples = nsamples;  %How Many Recorded Samples?
    mcmc_params.nthin = 1; % How Often is a Sample Recorded?
    mcmc_params.doparallel = runParallel; % Parallel Option
    mcmc_params.dic = 1;
end

% Ensure init0 is correct size
if ~isfield(mcmc_params, 'init0')
    for i=1:mcmc_params.nchains
        mcmc_params.init0(i) = struct;
    end
end
if ~exist('covariate', 'var')
    covariate = zeros(1, nSubj);
end


speTol = spe;
speTol(~spe) = 1e-5;
speTol(speTol==1) = 1-1e-5;
% if ~isempty(spe0)
%     spe0(~spe0) = 1e-5;
%     spe0(spe0==1) = 1-1e-5;
% end
% Assign variables to the observed nodes


model_file = ['fit_mbs_group_' modelNum '.txt'];


monitorparams = {'posneg_fact_fb_perc','posneg_fact_fb_mem',...
    'posneg_fact_conf_perc','posneg_fact_conf_mem','spe_est',...
    'v0_init', 'spe_global', 'dec_var', ...
    'beta_fb_lr', 'beta_conf_lr','beta_post_bias' };
datastruct = struct('corr', squeeze(corr), 'conf', squeeze(conf), 'feedback', ...
    squeeze(feedback), 'spe', speTol, ...
    'ntrials', ntrials, 'nblocks', nBlocks, 'nsubj', nSubj, 'task', task,...
    'cov', covariate);


% Use JAGS to Sample
try
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, model_file), ...
        mcmc_params.init0, ...
        'doparallel' , mcmc_params.doparallel, ...
        'nchains', mcmc_params.nchains,...
        'nburnin', mcmc_params.nburnin,...
        'nsamples', mcmc_params.nsamples, ...
        'thin', mcmc_params.nthin, ...
        'dic', mcmc_params.dic,...
        'monitorparams', monitorparams, ...
        'savejagsoutput' , 0 , ...
        'verbosity' , 1 , ...
        'cleanup' , 1 , ...
        'workingdir' , tmpfolder);
    toc
catch ME
    % Remove temporary directory if specified
    if exist('name','var')
        if exist(['../', tmpfolder],'dir')
            rmdir(['../', tmpfolder], 's');
        end
    end
    % Print the error message
    rethrow(ME);
end

% Remove temporary directory if specified
if newdircreated
    if exist(tmpfolder,'dir')
        %         rmdir(tmpfolder, 's');
    end
end

cd(cwd)
% stats.dic

fit.samples.beta_fb_lr = samples.beta_fb_lr;
fit.samples.beta_conf_lr = samples.beta_conf_lr;
if isfield(samples, 'beta_base')
    fit.samples.beta_base = samples.beta_base;
end
if isfield(samples, 'beta_post_bias')
    fit.samples.beta_post_bias = samples.beta_post_bias;
end
fit.samples.posneg_fact_fb_perc = samples.posneg_fact_fb_perc;
fit.samples.posneg_fact_fb_mem = samples.posneg_fact_fb_mem;
fit.samples.posneg_fact_conf_perc = samples.posneg_fact_conf_perc;
fit.samples.posneg_fact_conf_mem = samples.posneg_fact_conf_mem;

fit.means = stats.mean;

if doPlot
    figure,
    subplot(521), histogram((samples.posneg_fact_fb_perc))
    title('ΔLR feedback perc')
    subplot(522), histogram((samples.posneg_fact_conf_perc))
    title('ΔLR conf perc')
    subplot(523), histogram((samples.posneg_fact_fb_mem))
    title('ΔLR feedback mem')
    subplot(524), histogram((samples.posneg_fact_conf_mem))
    title('ΔLR conf mem')
    subplot(525), histogram((samples.beta_fb_lr(:,:,1)))
    title('beta(ΔLR) feedback')
    subplot(526), histogram((samples.beta_conf_lr(:,:,1)))
    title('beta(ΔLR) conf')
    subplot(528), histogram((samples.spe_global))
    title('spe mean prior')
    if isfield(samples, 'beta_base')
        subplot(527), histogram((samples.beta_base))
        title('beta base')
    end
    if isfield(samples, 'beta_post_bias')
        subplot(527), histogram((samples.beta_post_bias(:,:,1)))
        title('beta additive bias')
    end
    subplot(529), histogram((samples.v0_init(:,:,1)))
    title('spe variance prior perc')
    subplot(5,2,10), histogram((samples.v0_init(:,:,2)))
    title('spe variance prior mem')


end
fit.dic = stats.dic;
fit.spe_est = stats.mean.spe_est;

end
