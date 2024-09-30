function fit = fit_globalSPE(spe0,spe,corr,conf,feedback,task, ...
    ntrials,modelName,...
    nsamples, covariate, tmpfolder, fbBlock, doPlot)
%
%
%

sz = size(corr);
nSubj = sz(1);
% nTrials = sz(3);
nBlocks = sz(2);
% trialOrder = 1:nTrials;

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
    mcmc_params.doparallel = 1; % Parallel Option
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
if ~isempty(spe0)
    spe0(~spe0) = 1e-5;
    spe0(spe0==1) = 1-1e-5;
end
% Assign variables to the observed nodes


model_file = ['fit_mbs_group_' modelName '.txt'];

mean_conf = mean(conf,3, "omitnan");
mc1 = mean_conf(:,1:2).*(task(:,1:2)==0);
mc1(~mc1) = NaN;
mc2 = mean_conf(:,1:2).*(task(:,1:2)==1);
mc2(~mc2) = NaN;
meanconf(:,1) = mean(mc1,2,'omitnan');
meanconf(:,2) = mean(mc2,2,'omitnan');
meanconf(isnan(meanconf)) = 0;


monitorparams = {'posneg_fact_fb_perc','posneg_fact_fb_mem',...
    'posneg_fact_conf_perc','posneg_fact_conf_mem','spe_est',...
    'v0_init', 'spe_global', ....
    'beta_fb_lr', 'beta_conf_lr','beta_post_bias' };
monitorparams = {'beta_conf_lr'};
% datastruct = struct('corr', squeeze(corr), 'conf', squeeze(conf), 'feedback', ...
%     squeeze(feedback), 'spe', speTol, 'spe0', spe0,...
%     'ntrials', ntrials, 'nblocks', nBlocks, 'nsubj', nSubj, 'task', task,...
%     'cov', covariate, 'mean_conf', mean(mean(conf,2, "omitnan"),3, "omitnan"));
datastruct = struct(   'corr', squeeze(corr), 'conf', squeeze(conf), 'feedback', ...
    squeeze(feedback), 'spe', speTol, ...
    'ntrials', ntrials, 'nblocks', nBlocks, 'nsubj', nSubj, 'task', task,...
    'cov', covariate,  'meanconf', meanconf);

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
if isfield(samples, 'beta_localglobal')
    fit.samples.beta_localglobal = samples.beta_localglobal;
end
fit.samples.posneg_fact_fb_perc = samples.posneg_fact_fb_perc;
fit.samples.posneg_fact_fb_mem = samples.posneg_fact_fb_mem;
fit.samples.posneg_fact_conf_perc = samples.posneg_fact_conf_perc;
fit.samples.posneg_fact_conf_mem = samples.posneg_fact_conf_mem;
fit.samples.v0_init = samples.v0_init;
fit.samples.spe0    = samples.spe_global;

% fit.stats = stats;

if doPlot
    figure,
    subplot(521), histogram((samples.posneg_fact_fb_perc))
    title('del-LR fb perc')
    subplot(522), histogram((samples.posneg_fact_conf_perc))
    title('del-LR conf perc')
    subplot(523), histogram((samples.posneg_fact_fb_mem))
    title('del-LR fb mem')
    subplot(524), histogram((samples.posneg_fact_conf_mem))
    title('del-LR conf mem')
    subplot(525), histogram((samples.beta_fb_lr(:,:,1)))
    title('beta-fb-lr')
    subplot(526), histogram((samples.beta_conf_lr(:,:,1)))
    title('beta-conf-lr')
    subplot(527), histogram((samples.spe_global))
    title('spe global')
    if isfield(samples, 'beta_base')
        subplot(528), histogram((samples.beta_base))
        title('beta base')
    end
    if isfield(samples, 'beta_post_bias')
        subplot(528), histogram((samples.beta_post_bias(:,:,1)))
        title('beta post bias')
    end
    if isfield(samples, 'beta_localglobal')
        subplot(528), histogram((samples.beta_localglobal(:,:,1)))
        title('beta local-global')
    end
    subplot(529), histogram((samples.v0_init(:,:,1)))
    title('v0 perc')
    subplot(5,2,10), histogram((samples.v0_init(:,:,2)))
    title('v0 mem')


end
fit.dic = stats.dic;
fit.spe_est = stats.mean.spe_est;

end
