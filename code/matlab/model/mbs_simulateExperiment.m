
if ispc
    addpath('C:\Users\skatyal\OneDrive - University College London\Projects\Analysis_Such')
    addpath('C:\Users\skatyal\OneDrive - University College London\Projects\Experiments\meta_training_reanalysis\model_bias')
    addpath(genpath('C:\Users\skatyal\OneDrive - University College London\Projects\NeuroMatlabToolboxes/HMeta-d/'))

else
    addpath('/Users/skatyal/OneDrive - University College London/Projects/Experiments/meta_training_reanalysis/model_bias')
    addpath('/Users/skatyal/OneDrive - University College London/Projects/NeuroMatlabToolboxes/HMeta-d/Matlab/')
    addpath('/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/analysis/model/')
end

%% fit distribution to depression scores so we can randomly sample from them
% figure, histogram(mbsData.questAvg(:,1))

% phq = mbsData.questAvg(:,1);
% phq_dist= fitdist(phq,"Gamma");
% hold on
% histogram(random(phq_dist,300,1))
% % both halfnormal and gamma work well, though latter falls more sharply
% % similar to depression scores

% now make the distribution with the above parameter values, so the
% simulation doesn't require calling external functions
phq_dist= makedist("Gamma", "a",1.24972, "b",4.36897);
% hold on
% figure,histogram(random(phq_dist,300,1))

%% model metacognitive training v3 - all blocks

n_trials = 40; % number of trials in a session
n_subj = 600;

dprimes = [1.1, 0.05]; % two rows are the two tasks
mratios = [.8, .3; 1, .3];
c1 = 0;

spe_prior_mean = [.5 .01];
spe_prior_var = [.1 .05];

n_blocks = 6;

stairConvergence = .71;
p_feedback = [1/((1-stairConvergence)*n_trials), 9/((1-stairConvergence)*n_trials);...
    9/(stairConvergence*n_trials), 1/(stairConvergence*n_trials)]; % highest and lowest probability with which to give feedback on a trial type

subj_group = numel(1, n_subj);

groupFeedback = [...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0];
groupTask = [...
    0 0 0 0 0 0;...
    0 0 0 0 0 0;...
    0 1 0 1 0 1;...
    0 1 0 1 0 1;...
    1 1 1 1 1 1;...
    1 1 1 1 1 1;...
    0 1 1 0 1 0;...
    0 1 1 0 1 0];

fprintf('start... ')

stims = [-1 1]; % 2 stimuli

tmpfolder = 'tmpjags';
njagiters = 70;

mrdist_p = makedist('Normal','mu',mratios(1,1),'sigma',mratios(1,2));
mrdist_p = truncate(mrdist_p,0,2.5);
mrdist_m = makedist('Normal','mu',mratios(2,1),'sigma',mratios(2,2));
mrdist_m = truncate(mrdist_m,0,2.5);

dpdist = makedist('Normal','mu',dprimes(1),'sigma',dprimes(2));
dpdist = truncate(dpdist,0,2.5);

spe_mean_dist = makedist('Normal','mu',spe_prior_mean(1),'sigma',spe_prior_mean(2));
spe_mean_dist = truncate(spe_mean_dist,0,2.5);
spe_var_dist = makedist('Normal','mu',spe_prior_var(1),'sigma',spe_prior_var(2));
spe_var_dist = truncate(spe_var_dist,0,.3);

%%
modelNum = 6;
switch modelNum
    case 1
        % recover difference in learning rate
        beta_lr_fb = -.04:.02:.04;
        %         beta_lr_neg = -.05:.025:.05;
        beta_lr_conf = -.04:.02:.04;
        %         beta_lr_neg_conf = -.05:.025:.05;

        negpos_diff_fb = -.4:.2:.4;
        negpos_diff_conf = -.4:.2:.4;
        %         beta_lr_neg = -.05:.025:.05;

        beta_prior_mean = 0;
        beta_prior_var = 0;

        modelStr = num2str(5);

    case 2
        % does phq baseline corr with spe spuriously inject conf asymm.
        beta_lr_fb = 0;
        beta_lr_conf = 0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.2;

        beta_prior_mean = -.1;
        beta_prior_var = 0;

        modelStr = num2str(5);

    case 3
        % simulate model with only feedback asymmetry

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
        % simulate model with only confidence asymmetry

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
        % simulate model with noasymmetry

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = -.0;
        beta_lr_conf = -.0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = -.0;
        beta_prior_var = 0;
        beta_post_shift = -.0;

        %         modelStr = num2str(5);
    case 6
        % simulate model with additive asymmetry

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = 0;
        beta_lr_conf = 0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = 0;
        beta_prior_var = 0;
        beta_post_shift = -.005;

    case 7
        % simulate model with asymmetry in SPE prior

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = 0;
        beta_lr_conf = 0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = -.1;
        beta_prior_var = 0;
        beta_post_shift = 0;


    case 8
        % recover beta post bias

        beta_lr_fb = [-.02 .02];
        beta_lr_conf = [-.02 .02];
        beta_post_shift = [-.006 -.002:.001:.002 .006];

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;
        %         beta_lr_neg = -.05:.025:.05;

        beta_prior_mean = 0;
        beta_prior_var = 0;

        modelStr = num2str(7);

    case 9
        % simulate model with additive asymmetry with empirical values from
        % Exp 1

        n_subj = 5000;

        % simulate fitted data
        beta_lr_fb = 0;
        beta_lr_conf = 0;
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_prior_mean = 0;
        beta_prior_var = 0;
        beta_post_shift = 1.521635e-05;
    
    case 10
        % simulate model with additive asymmetry with empirical values from
        % Exp 1

        beta_lr_fb = 0;
        beta_lr_conf = 0;
        beta_post_shift = [-.005:.0025:.005];

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        beta_prior_mean = 0;
        beta_prior_var = 0;

        modelStr = num2str(10);
end
nlrfb = numel(beta_lr_fb);
nlr = numel(negpos_diff_fb);
nlrconf = numel(beta_lr_conf);
nlrpostbias = numel(beta_post_shift);
% nsp0mean = numel(beta_prior_mean);

%%

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

%%
initvals = [nlrf, nlrc, np, npc, npb];
for nlrf = initvals(1):nlr
    for nlrc = initvals(2) :nlr
        for np = initvals(3) :nlrfb
            for npc = initvals(4) :nlrconf
                for npb = initvals(5):nlrpostbias
                    %         for nsm = 1:nsp0mean
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
                            %                         spe(ns,nb) =
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
                            % in case of modelnum 3 just save the simulated data to
                            % make plots out of it in r
                            save(['data/exp1/simExp1_m' num2str(modelNum) '.mat'], ...
                                "spe", "phq", "task", "feedback", ...
                                "conf", "correct", "fbblock", 'group')

                        otherwise
                            fit = fit_mbsFeedback_group([], spe, correct, conf, feedback, ...
                                task, repmat(40, 1, 6), modelStr,njagiters, phq, tmpfolder);

                            fit_blr_fb(nlrf,nlrc,np,npc,npb) = mean(fit.samples.beta_fb_lr(:));
                            %         fit_lrneg(np,nn,npc,nnc) = mean(fit.samples.beta_lr_neg(:));
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

%%
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


    case 2
        figure, scatter(fit_betalr,beta_lr_fb)

    case {3, 4}
        figure, scatter(beta_prior_var, fit_betavar)

end

%%
%%
%%
%%
%% model metacognitive training v2 - all blocks

n_trials = 40; % number of trials in a session
n_subj = 600;

dprimes = [1.1, .1]; % two rows are the two tasks
mratios = [.8, .4; 1, .5];
c1 = 0;

spe_prior_mean = [.5 .1];
spe_prior_var = [.1 .1];

% lr_fb = [1;1]; %correct / incorrect

n_blocks = 6;

stairConvergence = .71;
p_feedback = [1/((1-stairConvergence)*n_trials), 9/((1-stairConvergence)*n_trials);...
    9/(stairConvergence*n_trials), 1/(stairConvergence*n_trials)]; % highest and lowest probability with which to give feedback on a trial type

% nlr = size(lr_fb,2);
subj_group = numel(1, n_subj);

groupFeedback = [...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0];
groupTask = [...
    0 0 0 0 0 0;...
    0 0 0 0 0 0;...
    0 1 0 1 0 1;...
    0 1 0 1 0 1;...
    1 1 1 1 1 1;...
    1 1 1 1 1 1;...
    0 1 1 0 1 0;...
    0 1 1 0 1 0];

fprintf('start... ')

stims = [-1 1]; % 2 stimuli

tmpfolder = 'tmpjags';
njagiters = 10;

mrdist_p = makedist('Normal','mu',mratios(1,1),'sigma',mratios(1,2));
mrdist_p = truncate(mrdist_p,0,2.5);
mrdist_m = makedist('Normal','mu',mratios(2,1),'sigma',mratios(2,2));
mrdist_m = truncate(mrdist_m,0,2.5);

dpdist = makedist('Normal','mu',dprimes(1),'sigma',dprimes(2));
dpdist = truncate(dpdist,0,2.5);

spe_mean_dist = makedist('Normal','mu',spe_prior_mean(1),'sigma',spe_prior_mean(2));
spe_mean_dist = truncate(spe_mean_dist,0,2.5);
spe_var_dist = makedist('Normal','mu',spe_prior_var(1),'sigma',spe_prior_var(2));
spe_var_dist = truncate(spe_var_dist,0,.3);

modelNum = 1;
switch modelNum
    case 1

        % recover difference in learning rate
        pos_diff_fb = 0;
        neg_diff_fb = 0;
        pos_diff_conf = .2 ;
        neg_diff_conf = .2;

        %         pos_diff_fb = -.4:.2:.4 ;
        %         neg_diff_fb = -.4:.2:.4 ;
        %         pos_diff_conf = -.4:.2:.4 ;
        %         neg_diff_conf = -.4:.2:.4 ;

        modelStr = num2str(4);


end
n_diff_fb = numel(pos_diff_fb);
n_diff_conf = numel(pos_diff_conf);

n_blr = numel(pos_diff_fb);
n_bvar = numel(pos_diff_conf);

%%
dprime = NaN(n_diff,n_diff_conf, n_subj, 2);
mratio = NaN(n_diff, n_diff_conf,n_subj, 2);

fit_pos_fb = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
fit_neg_fb = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
fit_pos_conf = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
fit_neg_conf = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
sim_pos_fb = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
sim_neg_fb = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
sim_pos_conf = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
sim_neg_conf = NaN(n_diff_fb, n_diff_conf,n_diff_fb, n_diff_conf,2);
model_iter = 0;

ntrials = repmat(40, 1, 6);

for nbv = 1:n_diff_fb
    for nblr = 1:n_diff_conf
        for nnp = 1:n_diff_fb
            for nnc = 1:n_diff_conf

                sim_pos_fb(nbv,nblr,nnp,nnc,1) = pos_diff_fb(nnp);
                sim_pos_fb(nbv,nblr,nnp,nnc,2) = pos_diff_fb(nnp);
                sim_neg_fb(nbv,nblr,nnp,nnc,1) = neg_diff_fb(nbv);
                sim_neg_fb(nbv,nblr,nnp,nnc,2) = neg_diff_fb(nbv);
                sim_pos_conf(nbv,nblr,nnp,nnc,1) = pos_diff_conf(nnc);
                sim_pos_conf(nbv,nblr,nnp,nnc,2) = pos_diff_conf(nnc);
                sim_neg_conf(nbv,nblr,nnp,nnc,1) = neg_diff_conf(nblr);
                sim_neg_conf(nbv,nblr,nnp,nnc,2) = neg_diff_conf(nblr);


                conf = NaN(n_subj, n_blocks, n_trials);
                feedback = conf;
                correct = conf;
                task = NaN(1, n_subj);
                spe = NaN(n_subj, n_blocks);
                perf = spe;
                spe0 = NaN(n_subj, n_blocks);
                spe0_var = task;
                phq = task;

                for ns = 1:n_subj

                    % assign one of 8 groups
                    subj_group(ns) = mod(ns,8)+1;

                    % task for this block
                    task = groupTask(subj_group(ns),:);

                    % sample mratio, dprime, spe prior mean and spe prior variance from
                    % a noisy distribution
                    mratio= [random(mrdist_p), random(mrdist_m)];
                    dprime= random(dpdist,1,2);
                    %                     mratio(nnp,nnc,ns,:) = [random(mrdist_p), random(mrdist_m)];
                    %                     dprime(nnp,nnc,ns,:) = random(dpdist,1,2);
                    phq(ns) = round(random(phq_dist));
                    spe0(ns,1) = random(spe_mean_dist);
                    spe0_var(ns) = random(spe_var_dist) + beta_prior_var(nbv)*phq(ns);

                    delta_lr_pos = 1 + (pos_diff_fb(nnp) + phq(ns)*beta_lr_fb(nblr))/2;
                    delta_lr_neg = 1 - (neg_diff_fb(nnp) + phq(ns)*beta_lr_fb(nblr))/2;
                    delta_lr_pos_conf = 1 + (pos_diff_conf(nnc) + phq(ns)*beta_lr_fb(nblr))/2;
                    delta_lr_neg_conf = 1 - (neg_diff_conf(nnc) + phq(ns)*beta_lr_fb(nblr))/2;

                    for nb = 1:n_blocks

                        metad = dprime(task(nb)+1)*mratio(task(nb)+1);

                        % current self-performance estimate
                        [aa,bb] = beta_mv2ab(...
                            spe0(ns,nb), spe0_var(ns));

                        for nt = 1:n_trials

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
                                    aa = aa + conf(ns,nb,nt) * delta_lr_pos ;
                                    bb = bb + (1-conf(ns,nb,nt))* delta_lr_neg;
                                end
                            else
                                % if no feedback, update with
                                % confidence
                                feedback(ns,nb,nt) = 0;
                                aa = aa + conf(ns,nb,nt) * delta_lr_pos_conf ;
                                bb = bb + (1-conf(ns,nb,nt))* delta_lr_neg_conf;
                            end
                        end
                        spe(ns,nb) = beta_ab2mv(aa,bb);
                        spe0(ns,nb+1) = spe(ns,nb);

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

                fit = fit_mbsFeedback_group([], spe, correct, conf, feedback, ...
                    task, ntrials, modelStr,njagiters, phq, tmpfolder);

                %                 figure,
                %                 subplot(121)
                %                 histogram(fit.samples.)
                fit_posneg_fb(nbv,nblr,nnp,nnc,1) = mean(fit.samples.posneg_fact_fb_perc(:));
                fit_posneg_fb(nbv,nblr,nnp,nnc,2) = mean(fit.samples.posneg_fact_fb_mem(:));
                fit_posneg_conf(nbv,nblr,nnp,nnc,1) = mean(fit.samples.posneg_fact_conf_perc(:));
                fit_posneg_conf(nbv,nblr,nnp,nnc,2) = mean(fit.samples.posneg_fact_conf_mem(:));
                fit_pos_fb(nbv,nblr,nnp,nnc,1) = pos_diff_fb(nnp);
                fit_pos_fb(nbv,nblr,nnp,nnc,2) = pos_diff_fb(nnp);
                fit_neg_fb(nbv,nblr,nnp,nnc,1) = neg_diff_fb(nnp);
                fit_neg_fb(nbv,nblr,nnp,nnc,2) = neg_diff_fb(nnp);
                fit_pos_conf(nbv,nblr,nnp,nnc,1) = pos_diff_conf(nnc);
                fit_pos_conf(nbv,nblr,nnp,nnc,2) = pos_diff_conf(nnc);
                fit_neg_conf(nbv,nblr,nnp,nnc,1) = neg_diff_conf(nnc);
                fit_neg_conf(nbv,nblr,nnp,nnc,2) = neg_diff_conf(nnc);

                %             fit_betalr(nbv,nblr,nnp,nnc) = mean(fit.samples.beta_lr(:));
                %             fit_betavar(nbv,nblr,nnp,nnc) = mean(fit.samples.beta_base_var(:));
                %             fit_betamean(nbv,nblr,nnp,nnc) = mean(fit.samples.beta_base(:));
            end
        end
        if n_diff_fb>1
        end
    end
end

%%
%% model metacognitive training v2 - if spe corr withAD spuriously injects asymm in conf

n_trials = 40; % number of trials in a session
n_subj = 600;

dprimes = [1.1, .1]; % two rows are the two tasks
mratios = [.8, .4; 1, .5];
c1 = 0;

spe_prior_mean = [.5 .1];
spe_prior_var = [.1 .1];

% lr_fb = [1;1]; %correct / incorrect

n_blocks = 6;

stairConvergence = .71;
p_feedback = [1/((1-stairConvergence)*n_trials), 9/((1-stairConvergence)*n_trials);...
    9/(stairConvergence*n_trials), 1/(stairConvergence*n_trials)]; % highest and lowest probability with which to give feedback on a trial type

% nlr = size(lr_fb,2);
subj_group = numel(1, n_subj);

groupFeedback = [...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0;...
    0 0 1 0 2 0;...
    0 0 2 0 1 0];
groupTask = [...
    0 0 0 0 0 0;...
    0 0 0 0 0 0;...
    0 1 0 1 0 1;...
    0 1 0 1 0 1;...
    1 1 1 1 1 1;...
    1 1 1 1 1 1;...
    0 1 1 0 1 0;...
    0 1 1 0 1 0];

fprintf('start... ')

stims = [-1 1]; % 2 stimuli

tmpfolder = 'tmpjags';
njagiters = 30;

mrdist_p = makedist('Normal','mu',mratios(1,1),'sigma',mratios(1,2));
mrdist_p = truncate(mrdist_p,0,2.5);
mrdist_m = makedist('Normal','mu',mratios(2,1),'sigma',mratios(2,2));
mrdist_m = truncate(mrdist_m,0,2.5);

dpdist = makedist('Normal','mu',dprimes(1),'sigma',dprimes(2));
dpdist = truncate(dpdist,0,2.5);

spe_mean_dist = makedist('Normal','mu',spe_prior_mean(1),'sigma',spe_prior_mean(2));
spe_mean_dist = truncate(spe_mean_dist,0,2.5);
spe_var_dist = makedist('Normal','mu',spe_prior_var(1),'sigma',spe_prior_var(2));
spe_var_dist = truncate(spe_var_dist,0,.3);

modelNum = 1;
switch modelNum
    case 1
        % recover difference in learning rate
        negpos_diff_fb = 0;
        negpos_diff_conf = 0;
        beta_lr_fb = 0;
        beta_prior_mean = 0.2;
        modelStr = '4';

end
n_diff_fb = numel(negpos_diff_fb);
n_diff_conf = numel(negpos_diff_fb);

n_blr = numel(beta_lr_fb);
n_bvar = numel(beta_prior_mean);

%%
dprime = NaN(n_diff_fb,n_diff_conf, n_subj, 2);
mratio = NaN(n_diff_fb, n_diff_conf,n_subj, 2);

fit_posneg_fb = NaN(n_bvar, n_blr, n_diff_fb, n_diff_conf,2);
fit_posneg_conf = NaN(n_bvar, n_blr, n_diff_fb, n_diff_conf,2);
sim_posneg_fb = NaN(n_bvar, n_blr, n_diff_fb, n_diff_conf,2);
sim_posneg_conf = NaN(n_bvar, n_blr, n_diff_fb, n_diff_conf,2);
% fit_betalr = NaN(n_bvar, n_blr, n_diff, n_diff_conf);
% fit_betavar = NaN(n_bvar, n_blr, n_diff, n_diff_conf);
% fit_betamean = NaN(n_bvar, n_blr, n_diff, n_diff_conf);
model_iter = 0;

ntrials = repmat(40, 1, 6);

for nbv = 1:n_bvar
    for nblr = 1:n_blr
        for nnp = 1:n_diff_fb
            for nnc = 1:n_diff_conf

                conf = NaN(n_subj, n_blocks, n_trials);
                feedback = conf;
                correct = conf;
                task = NaN(1, n_subj);
                spe = NaN(n_subj, n_blocks);
                perf = spe;
                spe0 = NaN(n_subj, n_blocks);
                spe0_var = task;
                phq = task;

                for ns = 1:n_subj

                    % assign one of 8 groups
                    subj_group(ns) = mod(ns,8)+1;

                    % task for this block
                    task = groupTask(subj_group(ns),:);

                    % sample mratio, dprime, spe prior mean and spe prior variance from
                    % a noisy distribution
                    mratio(nnp,nnc,ns,:) = [random(mrdist_p), random(mrdist_m)];
                    dprime(nnp,nnc,ns,:) = random(dpdist,1,2);
                    phq(ns) = round(random(phq_dist));
                    spe0(ns,1) = random(spe_mean_dist)+ beta_prior_mean(nbv)*phq(ns);
                    spe0_var(ns) = random(spe_var_dist) ;

                    sim_posneg_fb(nbv,nblr,nnp,nnc,1) = negpos_diff_fb(nnp);
                    sim_posneg_fb(nbv,nblr,nnp,nnc,2) = negpos_diff_fb(nnp);
                    sim_posneg_conf(nbv,nblr,nnp,nnc,1) = negpos_diff_conf(nnc);
                    sim_posneg_conf(nbv,nblr,nnp,nnc,2) = negpos_diff_conf(nnc);

                    delta_lr_pos = 1 + (negpos_diff_fb(nnp) + phq(ns)*beta_lr_fb(nblr))/2;
                    delta_lr_neg = 1 - (negpos_diff_fb(nnp) + phq(ns)*beta_lr_fb(nblr))/2;
                    delta_lr_pos_conf = 1 + (negpos_diff_conf(nnc) + phq(ns)*beta_lr_fb(nblr))/2;
                    delta_lr_neg_conf = 1 - (negpos_diff_conf(nnc) + phq(ns)*beta_lr_fb(nblr))/2;

                    for nb = 1:n_blocks

                        metad = dprime(nnp,ns,task(nb)+1)*mratio(nnp,ns,task(nb)+1);

                        % current self-performance estimate
                        [aa,bb] = beta_mv2ab(...
                            spe0(ns,nb), spe0_var(ns));

                        for nt = 1:n_trials

                            stim = stims((rand>.5)+1); % choose 1 of 2 stimuli

                            X = normrnd(stim * dprime(nnp,ns,task(nb)+1)/2, 1) ;

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
                                    aa = aa + conf(ns,nb,nt) * delta_lr_pos ;
                                    bb = bb + (1-conf(ns,nb,nt))* delta_lr_neg;
                                end
                            else
                                % if no feedback, update with
                                % confidence
                                feedback(ns,nb,nt) = 0;
                                aa = aa + conf(ns,nb,nt) * delta_lr_pos_conf ;
                                bb = bb + (1-conf(ns,nb,nt))* delta_lr_neg_conf;
                            end
                        end
                        spe(ns,nb) = beta_ab2mv(aa,bb);
                        spe0(ns,nb+1) = spe(ns,nb);

                    end

                    perf = mean(correct,3);

                end
                include_subj = all(perf>.6 & perf<.85,2) & all(spe<1 & spe>0,2);
                sum(include_subj)

                model_iter = model_iter+1;
                fprintf('model iteration: %g\n', model_iter)

                spe = spe(include_subj,:);
                correct = correct(include_subj,:,:);
                conf = conf(include_subj,:,:);
                feedback = feedback(include_subj,:,:);
                task = groupTask(subj_group(include_subj),:);
                phq = phq(include_subj);

                fit = fit_mbsFeedback_group([], spe, correct, conf, feedback, ...
                    task, ntrials, modelStr,njagiters, phq, tmpfolder);

                fit_posneg_fb(nbv,nblr,nnp,nnc,1) = mean(fit.samples.posneg_fact_fb_perc(:));
                fit_posneg_fb(nbv,nblr,nnp,nnc,2) = mean(fit.samples.posneg_fact_fb_mem(:));
                fit_posneg_conf(nbv,nblr,nnp,nnc,1) = mean(fit.samples.posneg_fact_conf_perc(:));
                fit_posneg_conf(nbv,nblr,nnp,nnc,2) = mean(fit.samples.posneg_fact_conf_mem(:));

                %             fit_betalr(nbv,nblr,nnp,nnc) = mean(fit.samples.beta_lr(:));
                %             fit_betavar(nbv,nblr,nnp,nnc) = mean(fit.samples.beta_base_var(:));
                %             fit_betamean(nbv,nblr,nnp,nnc) = mean(fit.samples.beta_base(:));
            end
        end
        if n_diff_fb>1
        end
    end
end
