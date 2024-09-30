
addpath('~/OneDrive - University of Copenhagen/Projects/Experiments/meta_training_reanalysis/model_bias')
addpath('~/OneDrive - University of Copenhagen/Projects/NeuroMatlabToolboxes/HMeta-d/Matlab/')
addpath('~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/analysis/model/')

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
phq_dist= makedist("Gamma", "a",1.3, "b",4.4);
hold on
figure,histogram(random(phq_dist,1000,1))

%% model metacognitive training v3 - all blocks

n_trials = 40; % number of trials in a session
n_subj = 600;
% n_subj = 2000;

dprimes = [1.1, 0.05]; % two rows are the two tasks
mratios = [.8, .3; 1, .3];
% mratios = [.5, .5; .8, .5];
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

tmpfolder = 'tmpjags5';
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

modelNum = 104;


switch modelNum
    case 1
        % recover difference in learning rates for the feedback and
        % confidence asymmetry regression slope parameters along with
        % regression  intercepts
        beta_lr_fb = -.04:.02:.04;
        beta_lr_conf = -.04:.02:.04;
        beta_post_shift = 0;

        negpos_diff_fb = -.4:.2:.4;
        negpos_diff_conf = -.4:.2:.4;

        modelStr = num2str(8);

    case 2
        % recover the additive bias regression slope parameter along with
        % the regression intercept
        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        beta_lr_fb = 0;
        beta_lr_conf = 0;
        beta_post_shift = -.00005:.000025:.00005;

        modelStr = num2str(7);


    case 3
        % simulate model with only feedback asymmetry

        % n_subj = 1000;

        % simulate fitted data
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_lr_fb = -.1;
        beta_lr_conf = 0;
        beta_post_shift = 0;

        % beta_prior_mean = -.0;
        % beta_prior_var = 0;


    case 4
        % simulate model with only confidence asymmetry

        % n_subj = 1000;

        % simulate fitted data
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_lr_fb = 0;
        beta_lr_conf = -.08;
        beta_post_shift = 0;


    case 5
        % simulate model with noasymmetry

        % n_subj = 1000;

        % simulate fitted data
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_lr_fb = -.0;
        beta_lr_conf = -.0;
        beta_post_shift = -.0;

        % beta_prior_mean = -.0;
        % beta_prior_var = 0;

        %         modelStr = num2str(5);
    case 6
        % simulate model with additive asymmetry

        % n_subj = 1000;

        % simulate fitted data
        negpos_diff_fb = -.0;
        negpos_diff_conf = -.1;
        % negpos_diff_fb = -.3;
        % negpos_diff_conf = .3;

        beta_lr_fb = 0;
        beta_lr_conf = 0;
        beta_post_shift = -.01;

    case 7
        % simulate model with both feedback and confidence asymmetry
        % regression slopes

        % n_subj = 5000;

        % simulate fitted data
        negpos_diff_fb = 0;
        negpos_diff_conf = 0;

        beta_lr_fb = -.1;
        beta_lr_conf = -.3;
        beta_post_shift = 0;


    case 9
        % simulate model with additive asymmetry with empirical values from
        % Exp 1

        % n_subj = 5000;

        % simulate fitted data
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_lr_fb = 0;
        beta_lr_conf = 0;
        beta_post_shift = 1.521635e-05;


    case 10
        % recover all 3 betas

        beta_lr_fb = [-.02 0 .02];
        beta_lr_conf = [-.02 0 .02];
        beta_post_shift = [-.005 0 .005];

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        modelStr = num2str(12);

    case 11
        % simulate conf distortion recover conf+respbias

        beta_lr_fb = 0;
        beta_lr_conf = [-.02 0 .02];
        beta_post_shift = 0;

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        modelStr = num2str(10);

    case 12
        % simulate conf distortion recover conf+respbias

        beta_lr_fb = 0;
        beta_lr_conf = [-.02 0 .02];
        beta_post_shift = [-.02 0 .02];

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        modelStr = num2str(10);

    case 13
        % simulate fb distortion recover respbias

        beta_lr_fb = [-.02 0 .02];
        beta_lr_conf = 0;
        beta_post_shift = 0;

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        modelStr = num2str(7);

    case 14
        % simulate no distortion recover respbias

        beta_lr_fb = 0;
        beta_lr_conf = 0;
        beta_post_shift = 0;

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        modelStr = num2str(7);

    case 15
        % simulate conf distortion recover respbias

        beta_lr_fb = 0;
        beta_lr_conf = 0;
        beta_post_shift = 0;

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        modelStr = num2str(7);

    case 16
        % simulate respbias distortion recover conf

        beta_lr_fb = 0;
        beta_lr_conf = 0;
        beta_post_shift = [-.0005 0 .0005];

        negpos_diff_fb = -.25:.25:.25;
        negpos_diff_conf = -.25:.25:.25;

        modelStr = num2str(6);

    case 17
        % simulate model with additive asymmetry with empirical values from
        % Exp 1

        % n_subj = 5000;

        % simulate fitted data
        negpos_diff_fb = 0;
        negpos_diff_conf = -.1;

        beta_lr_fb = -.01;
        beta_lr_conf = 0;
        beta_post_shift = 0;


    case 104
        % simulate model with only confidence asymmetry

        % n_subj = 1000;

        % simulate fitted data
        negpos_diff_fb = .0;
        negpos_diff_conf = 0;

        beta_lr_fb = 0;
        beta_lr_conf = -.08;
        beta_post_shift = 0;

        modelStr = num2str(4);


end
nlrfb = numel(beta_lr_fb);
nlr = numel(negpos_diff_fb);
nlrconf = numel(beta_lr_conf);
nlrpostbias = numel(beta_post_shift);

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
                        spe0(ns,1) = random(spe_mean_dist);% + beta_prior_mean*phq(ns);
                        spe0_var(ns) = random(spe_var_dist);
                        % confbias(ns) = random(conf_bias_dist);
                        % confbias(ns) = 0;

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

                                % normalise conf to 0-1
                                conf(ns, nb, nt) = 2*conf(ns, nb, nt) - 1;
                                % if conf(ns, nb, nt) > 1, conf(ns, nb, nt) = 1; end
                                % if conf(ns, nb, nt) < 0, conf(ns, nb, nt) = 0; end
                                
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
                            if aa<=0, aa=1e-5; end 
                            if bb<=0, bb=1e-5; end
                            spe0(ns,nb+1) = beta_ab2mv(aa,bb);
                            spe(ns,nb)= spe0(ns,nb+1) + beta_post_shift(npb)*phq(ns); % with no learning to the next block
                        end

                        perf = mean(correct,3);

                    end
                    include_subj = all(perf>.6 & perf<=.85,2) & all(spe<1 & spe>0,2);
                    % include_subj = true(1, ns);
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
                        case {3,4,5,6,7,9,17,104}
                            % in case of modelnum 3-9 just save the simulated data to
                            % make plots out of it in r
                            save(['data/simulation/simExp_m' num2str(modelNum) '.mat'], ...
                                "spe", "phq", "task", "feedback", ...
                                "conf", "correct", "fbblock", 'group')
                            fprintf('saved %s\n', ['data/simulation/simExp_m' num2str(modelNum) '.mat'])
                        otherwise
                            fprintf('fb: %g,  conf: %g, betafb: %g, betaconf: %g, betaadd: %g\n', ...
                                nlrf, nlrc, np, npc, npb)

                            fit = fit_globalSPE2([], spe, correct, conf, feedback, ...
                                task, repmat(40, 1, 6), modelStr, njagiters, phq, tmpfolder);

                            fit_lr_fb(nlrf,nlrc,np,npc,npb) = mean(fit.samples.posneg_fact_fb_perc(:));
                            fit_lr_conf(nlrf,nlrc,np,npc,npb) = mean(fit.samples.posneg_fact_conf_perc(:));

                            fit_blr_fb(nlrf,nlrc,np,npc,npb) = mean(fit.samples.beta_fb_lr(:));
                            fit_bpostbias(nlrf,nlrc,np,npc,npb) = mean(fit.samples.beta_post_bias(:));
                            fit_blr_conf(nlrf,nlrc,np,npc,npb) = mean(fit.samples.beta_conf_lr(:));
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
        save('data/paramRecovery_fbconf.mat')
    case 2
        save('data/paramRecovery_postbias_smallvals.mat')
    case 10
        save('data/paramRecovery_fbconfpost.mat')
    case 11
        save('data/paramRecovery_simConf_recConfplusBias.mat')
    case 12
        save('data/paramRecovery_simConfplusBias_recConfplusBias.mat')
    case 13
        save('data/paramRecovery_simFb_recBias.mat')
    case 14
        save('data/paramRecovery_simNoDist_recBias.mat')
    case 15
        save('data/paramRecovery_simConf_recBias.mat')
end

%%
figure, scatter( sim_blr_conf(:), fit_blr_conf(:))
figure, scatter( sim_blr_fb(:), fit_blr_fb(:))
figure, scatter( sim_bpostbias(:), fit_bpostbias(:))
line([-.00005 .00005], [-.00005, .00005], 'Color', 'k', 'LineWidth', 2)
xlabel('Simulated additive bias', 'FontSize', 16)
ylabel('Recovered additive bias', 'FontSize', 16)

figure, scatter( sim_blr_conf(:), fit_blr_fb(:))
figure, scatter( sim_blr_fb(:), fit_blr_conf(:))
figure, scatter( sim_blr_conf(:), fit_bpostbias(:))
figure, scatter( sim_blr_fb(:), fit_bpostbias(:))

%% this cell plots comparisons of 2 (fb+conf) vs. 3 (fb+conf+bias) regression parameter model recovery - used to respond to Reviewer 3

close all

load('data/paramRecovery_fbconf.mat')
figure(1), scatter( sim_blr_conf(abs(sim_blr_conf)<.03)-.001, ...
    fit_blr_conf(abs(sim_blr_conf)<.03))
figure(2), scatter( sim_blr_fb(abs(sim_blr_fb)<.03)-.001, ...
    fit_blr_fb(abs(sim_blr_fb)<.03))

load('data/paramRecovery_fbconfpost.mat')
figure(1), hold on, scatter( sim_blr_conf(abs(sim_blr_conf)<.03)+.001, ...
    fit_blr_conf(abs(sim_blr_conf)<.03))
line([-.03 .03], [-.03, .03], 'Color', 'k', 'LineWidth', 2)
xlabel('Simulated β-confidence', 'FontSize', 16)
ylabel('Recovered β-confidence', 'FontSize', 16)

figure(2), hold on, scatter( sim_blr_fb(abs(sim_blr_fb)<.03)+.001, ...
    fit_blr_fb(abs(sim_blr_fb)<.03))
line([-.03 .03], [-.03, .03], 'Color', 'k', 'LineWidth', 2)
xlabel('Simulated β-feedback', 'FontSize', 16)
ylabel('Recovered β-feedback', 'FontSize', 16)


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


end

