% 
% if ispc
%     addpath('C:\Users\skatyal\OneDrive - University College London\Projects\Analysis_Such')
%     addpath('C:\Users\skatyal\OneDrive - University College London\Projects\Experiments\meta_training_reanalysis\model_bias')
%     addpath(genpath('C:\Users\skatyal\OneDrive - University College London\Projects\NeuroMatlabToolboxes/HMeta-d/'))
% 
% else
%     addpath('/Users/skatyal/OneDrive - University College London/Projects/Experiments/meta_training_reanalysis/model_bias')
%     addpath('/Users/skatyal/OneDrive - University College London/Projects/NeuroMatlabToolboxes/HMeta-d/Matlab/')
%     addpath('/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/analysis/model/')
% end
%% Script for simulating model data for the study
%% Katyal, Huys, Dolan, Fleming - How underconfidence is maintained in anxiety and depression
%%
% BEFORE PROCEEDING: Set current directory to the base directory of
% downloaded repo from github.com/sucharitk/confidence-distortion-AD

%% 1: add paths

clear
close all
addpath(genpath('code'))


%% 2: initialise parameters

n_trials = 40; % number of trials in a block
n_subj = 600; % number of subjects (before applying exclusion criteria)

dprimes = [1.1, 0.05]; % population mean and variance: corresponds to ~71% performance
mratios = [.8, .3; 1, .3]; % population mean and variance. two rows are for the two tasks, mratio: memory > perception
c1 = 0;

spe_prior_mean = [.5 .01]; % values similar to observed values
spe_prior_var = [.1 .05];

n_blocks = 6;

% define a distribution of phq (depression) scores with parameters obtained
% by fitting a Gamma distribution to the observed values in Exp 1
phq_dist= makedist("Gamma", "a",1.24972, "b",4.36897);

% expected accuracy based on the staircase
stairConvergence = .71;
% probabilities with which to give feedback on correct and incorrect trials
% based on expected accuracy for positive and negative blocks
p_feedback = [1/((1-stairConvergence)*n_trials), 9/((1-stairConvergence)*n_trials);...
    9/(stairConvergence*n_trials), 1/(stairConvergence*n_trials)]; 

% there will be 8 groups of participants with different orders of feedback
% and tasks
subj_group = numel(1, n_subj);
groupFeedback = [...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0]; % 1-pos , 2-neg feedback
groupTask = [...
    0 0 0 0 0 0;...
    0 0 0 0 0 0;...
    0 1 0 1 0 1;...
    0 1 0 1 0 1;...
    1 1 1 1 1 1;...
    1 1 1 1 1 1;...
    0 1 1 0 1 0;...
    0 1 1 0 1 0]; % two tasks, 0-perc, 1-memory


stims = [-1 1]; % 2 choices

tmpfolder = 'tmpjags';
njagiters = 70;

% distributions of mratios
mrdist_p = makedist('Normal','mu',mratios(1,1),'sigma',mratios(1,2));
mrdist_p = truncate(mrdist_p,0,2.5);
mrdist_m = makedist('Normal','mu',mratios(2,1),'sigma',mratios(2,2));
mrdist_m = truncate(mrdist_m,0,2.5);

% distributions of dprimes
dpdist = makedist('Normal','mu',dprimes(1),'sigma',dprimes(2));
dpdist = truncate(dpdist,0,2.5);

% distributions of spe mean and variance priors
spe_mean_dist = makedist('Normal','mu',spe_prior_mean(1),'sigma',spe_prior_mean(2));
spe_mean_dist = truncate(spe_mean_dist,0,2.5);
spe_var_dist = makedist('Normal','mu',spe_prior_var(1),'sigma',spe_prior_var(2));
spe_var_dist = truncate(spe_var_dist,0,.3);

fprintf('start... ')

%% 3: various modes of running the model

modelNum = 6;

switch modelNum

    %%% following cases used for plotting model recovery
    case 1
        % recover ∆LR and Beta(∆LR) for feedback and confidence
        negpos_diff_fb = -.4:.2:.4; % ∆LR used based on observed values
        negpos_diff_conf = -.4:.2:.4;

        beta_lr_fb = -.04:.02:.04;
        beta_lr_conf = -.04:.02:.04;

        beta_prior_mean = 0;
        beta_prior_var = 0;

        modelStr = num2str(8);

    case 2
        % recover Beta for additive bias
        negpos_diff_fb = -.25:.25:.25; % ∆LR used based on observed values
        negpos_diff_conf = -.25:.25:.25;

        beta_lr_fb = [-.02 .02];
        beta_lr_conf = [-.02 .02];
        beta_post_shift = [-.006 -.002:.001:.002 .006];

        beta_prior_mean = 0;
        beta_prior_var = 0;

        modelStr = num2str(7);

        %%% following cases used for plotting model predictions
    case 3
        % simulate model with only feedback distortion

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = -.2;
        beta_lr_conf = -.0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = -.0;
        beta_prior_var = 0;
        beta_post_shift = -.0;


    case 4
        % simulate model with only confidence distortion

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = -.0;
        beta_lr_conf = -.04;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = -.0;
        beta_prior_var = 0;
        beta_post_shift = -.0;


    case 5
        % simulate model with no distortion

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = -.0;
        beta_lr_conf = -.0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = -.0;
        beta_prior_var = 0;
        beta_post_shift = -.0;

    case 6
        % simulate model with additive bias

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = 0;
        beta_lr_conf = 0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = 0;
        beta_prior_var = 0;
        beta_post_shift = -.005;


    case 9
        % simulate model with additive asymmetry with empirical values from
        % Exp 1 - used in Supp Figure 6

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = 0;
        beta_lr_conf = 0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = 0;
        beta_prior_var = 0;
        beta_post_shift = 1.521635e-05;
    
end
nlrfb = numel(beta_lr_fb);
nlr = numel(negpos_diff_fb);
nlrconf = numel(beta_lr_conf);
nlrpostbias = numel(beta_post_shift);
% nsp0mean = numel(beta_prior_mean);

%% 4: initialise variables

fit_lr_fb = NaN(nlr, nlr, nlrfb, nlrconf, nlrpostbias);
fit_lr_conf = fit_lr_fb;
fit_blr_fb = fit_lr_fb;
fit_blr_conf = fit_lr_fb;
fit_bpostbias = fit_lr_fb;
fit_bpostbias_median = fit_lr_fb;
sim_lr_fb = fit_lr_fb;
sim_lr_conf = fit_lr_conf;
sim_blr_fb = fit_lr_fb;
sim_blr_conf = fit_lr_fb;
sim_bpostbias = fit_lr_fb;

nlrf = 1; nlrc = 1; np = 1; npc = 1; npb = 1;
model_iter = 0;

%% 5: run the simulation and recovery

initvals = [nlrf, nlrc, np, npc, npb];
for nlrf = initvals(1):nlr
    for nlrc = initvals(2) :nlr
        for np = initvals(3) :nlrfb
            for npc = initvals(4) :nlrconf
                for npb = initvals(5):nlrpostbias

                    conf = NaN(n_subj, n_blocks, n_trials);
                    feedback = conf;
                    correct = conf;

                    spe = NaN(n_subj, n_blocks);
                    perf = spe;
                    spe0 = NaN(n_subj, n_blocks);
                    spe0_var = NaN(1,n_subj);
                    phq = spe0_var;

                    for ns = 1:n_subj

                        % assign one of 8 groups
                        subj_group(ns) = mod(ns,8)+1;

                        % task for this block
                        task = groupTask(subj_group(ns),:);

                        % sample mratio, dprime, spe prior mean and spe prior variance from
                        % a noisy distribution
                        mratio = [random(mrdist_p), random(mrdist_m)];
                        dprime = random(dpdist,1,2);
                        phq(ns) = round(random(phq_dist));
                        spe0(ns,1) = random(spe_mean_dist) + beta_prior_mean*phq(ns);
                        spe0_var(ns) = random(spe_var_dist);

                        sim_blr_fb(nlrf,nlrc,np,npc,npb) = beta_lr_fb(np);
                        sim_blr_conf(nlrf,nlrc,np,npc,npb) = beta_lr_conf(npc);
                        sim_lr_fb(nlrf,nlrc,np,npc,npb) = negpos_diff_fb(nlrf);
                        sim_lr_conf(nlrf,nlrc,np,npc,npb) = negpos_diff_conf(nlrc);
                        sim_bpostbias(nlrf,nlrc,np,npc,npb) = beta_post_shift(npb);

                        delta_lr_pos = 1 + (negpos_diff_fb(nlrf) + phq(ns)*beta_lr_fb(np))/2;
                        delta_lr_neg = 1 - (negpos_diff_fb(nlrf) + phq(ns)*beta_lr_fb(np))/2;
                        delta_lr_pos_conf = 1 + (negpos_diff_conf(nlrc) + phq(ns)*beta_lr_conf(npc))/2;
                        delta_lr_neg_conf = 1 - (negpos_diff_conf(nlrc) + phq(ns)*beta_lr_conf(npc))/2;

                        for nb = 1:n_blocks

                            metad = dprime(task(nb)+1)*mratio(task(nb)+1);

                            % current self-performance estimate
                            [aa,bb] = beta_mv2ab(...
                                spe0(ns,nb), spe0_var(ns));

                            for nt = 1:n_trials
                                % simulate individual trials from an SDT
                                % model

                                stim = stims((rand>.5)+1); % choose 1 of 2 stimuli

                                X = normrnd(stim * dprime(task(nb)+1)/2, 1) ;

                                Xoff = X - c1;

                                if Xoff < 0
                                    resp = -1; % choose left
                                else
                                    resp = 1; % choose right
                                end

                                correct(ns, nb, nt) = resp == stim;

                                % for simulating meta-level sensitivity generate a second sample
                                % from meta-d' which is in the same direction from c1 as X
                                X2 = normrnd(stim * metad/2, 1) ;
                                while sign(X2-c1) ~= sign(Xoff)
                                    X2 = normrnd(stim * metad/2, 1) ;
                                end
                                Xoff = X2 - c1;

                                % compute confidence as probability of
                                % being correct
                                if resp > 0 % chose right
                                    conf(ns, nb, nt) = normpdf(Xoff, metad/2, 1) / ...
                                        (normpdf(Xoff, metad/2,1) + ...
                                        normpdf(Xoff, -metad/2,1)) ;
                                else               % chose a = -1 left
                                    conf(ns, nb, nt) = normpdf(Xoff, -metad/2, 1) / ...
                                        (normpdf(Xoff, metad/2, 1) + ...
                                        normpdf(Xoff, -metad/2, 1)) ;
                                end


                                %%% feedback module
                                if groupFeedback(subj_group(ns), nb)
                                    % if this block involves feedback
                                    feedback(ns,nb,nt) = (rand < ...
                                        p_feedback(correct(ns,nb,nt)+1, ...
                                        groupFeedback(subj_group(ns), nb)))...
                                        * (correct(ns,nb,nt)*2-1);

                                    if feedback(ns,nb,nt)
                                        if correct(ns,nb,nt)
                                            aa = aa + delta_lr_pos;
                                        else
                                            bb = bb + delta_lr_neg;
                                        end
                                    else
                                        aa = aa + conf(ns,nb,nt) * delta_lr_pos_conf ;
                                        bb = bb + (1-conf(ns,nb,nt))* delta_lr_neg_conf;
                                    end
                                else
                                    % if no feedback, update with
                                    % confidence
                                    feedback(ns,nb,nt) = 0;
                                    aa = aa + conf(ns,nb,nt) * delta_lr_pos_conf ;
                                    bb = bb + (1-conf(ns,nb,nt))* delta_lr_neg_conf;
                                end
                            end
                            spe0(ns,nb+1) = beta_ab2mv(aa,bb);
                            spe(ns,nb)= spe0(ns,nb+1) + beta_post_shift(npb)*phq(ns); % with no learning to the next block
                        end

                        perf = mean(correct,3);

                    end
                    include_subj = all(perf>.6 & perf<=.85,2) & all(spe<1 & spe>0,2);
                    sum(include_subj)

                    model_iter = model_iter+1;
                    fprintf('model iteration: %g\n', model_iter)

                    spe = spe(include_subj,:);
                    correct = correct(include_subj,:,:);
                    conf = conf(include_subj,:,:);
                    feedback = feedback(include_subj,:,:);
                    task = groupTask(subj_group(include_subj),:);
                    phq = phq(include_subj);
                    fbblock = groupFeedback(subj_group(include_subj),:);
                    group = subj_group(include_subj);

                    switch modelNum
                        case {3,4,5,6,7,9}
                            % just save the simulated data to
                            % make plots out of it in r (no recovery)
                            save(['data/exp1/simExp1_m' num2str(modelNum) '.mat'], ...
                                "spe", "phq", "task", "feedback", ...
                                "conf", "correct", "fbblock", 'group')

                        otherwise
                            % fit the model
                            fit = fit_mbsFeedback_group([], spe, correct, conf, feedback, ...
                                task, repmat(40, 1, 6), modelStr,njagiters, phq, tmpfolder);

                            fit_blr_fb(nlrf,nlrc,np,npc,npb) = mean(fit.samples.beta_fb_lr(:));
                            fit_blr_conf(nlrf,nlrc,np,npc,npb) = mean(fit.samples.beta_conf_lr(:));
                            fit_lr_fb(nlrf,nlrc,np,npc,npb) = mean(fit.samples.posneg_fact_fb_perc(:));
                            fit_lr_conf(nlrf,nlrc,np,npc,npb) = mean(fit.samples.posneg_fact_conf_perc(:));
                            fit_bpostbias(nlrf,nlrc,np,npc,npb) = mean(fit.samples.beta_post_bias(:));
                    end
                end
                initvals(5) = 1;
            end
            initvals(4) = 1;
        end
        initvals(3) = 1;
    end
    initvals(2) = 1;
end

%% 6: plot the models
% only implemented here for the first case, other cases are plotted using
% the R code mbsAnalyses_model.R
switch modelNum
    case 1
        fontSize = 15;
        figure,subplot(221), hold on
        line([min(sim_lr_fb(:)) max(sim_lr_fb(:))],...
            [min(sim_lr_fb(:)) max(sim_lr_fb(:))],'Color', 'k', 'LineWidth',3)
        scatter(sim_lr_fb(:), fit_blr_fb(:))
        ax = gca;
        ax.FontSize = fontSize;
        title('∆LR feedback','FontSize',16,'FontWeight','Normal')

        subplot(222), hold on,
        line([min(sim_lr_conf(:)) max(sim_lr_conf(:))],...
            [min(sim_lr_conf(:)) max(sim_lr_conf(:))],'Color','k', 'LineWidth',3)
        scatter(sim_lr_conf(:), fit_blr_conf(:))
        ax = gca;
        ax.FontSize = fontSize;
        title('∆LR confidence','FontSize',16,'FontWeight','Normal')
        subplot(223), hold on,
        line([min(sim_blr_fb(:)) max(sim_blr_fb(:))],...
            [min(sim_blr_fb(:)) max(sim_blr_fb(:))], 'Color','k','LineWidth',3)
        scatter(sim_blr_fb(:), fit_lr_fb(:))
        ax = gca;
        ax.FontSize = fontSize;
        title('Beta (∆LR feedback)','FontSize',16,'FontWeight','Normal')

        subplot(224), hold on,
        line([min(sim_blr_conf(:)) max(sim_blr_conf(:))],...
            [min(sim_blr_conf(:)) max(sim_blr_conf(:))], 'Color','k','LineWidth',3)
        scatter(sim_blr_conf(:), fit_lr_conf(:))
        ax = gca;
        ax.FontSize = fontSize;
        title('Beta (∆LR confidence)','FontSize',16,'FontWeight','Normal')



end


