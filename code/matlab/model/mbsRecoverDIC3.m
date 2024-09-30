
addpath('~/OneDrive - University of Copenhagen/Projects/Experiments/meta_training_reanalysis/model_bias')
addpath('~/OneDrive - University of Copenhagen/Projects/NeuroMatlabToolboxes/HMeta-d/Matlab/')
addpath('~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/analysis/model/')

%% fit distribution to depression scores so we can randomly sample from them

phq_dist= makedist("Gamma", "a",1.3, "b",4.4);
hold on
figure,histogram(random(phq_dist,1000,1))

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

tmpfolder = 'tmpjags3';
njagiters = 120;

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

negpos_diff_fb = -.4:.2:.4;
negpos_diff_conf = -.4:.2:.4;

% betas = [-.04 -.02 .02 .04];
% betas = [-.03 -.01];
betas = [-.04 -.03 -.02 -.01];

betatypes = eye(3);
nbetatypes = 3;
model_strs = [5,6,7];

nlr = numel(negpos_diff_fb);
nbeta = numel(betas);

%%

sim_lr_fb = NaN(nlr, nlr, nbeta, nbetatypes);
sim_lr_conf = sim_lr_fb;
sim_beta = sim_lr_fb;
sim_betasimtype = sim_lr_fb;
dic_vals = NaN(nlr, nlr, nbeta, nbetatypes, nbetatypes);

ilrf = 1; ilrc = 1; ib = 1; ibstsim = 1; ibstrec = 1;
% model_iter = 0;

%%
initvals = [ibstsim, ib, ilrf, ilrc];
% for ibstsim = initvals(1) :nbetatypes
ibstsim = 3;
    for ib = initvals(2) :nbeta
        for ilrf = initvals(3):nlr
            for ilrc = initvals(4):nlr

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

                    sim_lr_fb(ilrf,ilrc,ib,ibstsim) = negpos_diff_fb(ilrf);
                    sim_lr_conf(ilrf,ilrc,ib,ibstsim) = negpos_diff_conf(ilrc);
                    sim_beta(ilrf,ilrc,ib,ibstsim) = betas(ib);
                    sim_betasimtype(ilrf,ilrc,ib,ibstsim) = ib;

                    delta_lr_pos = 1 + (negpos_diff_fb(ilrf) + ...
                        phq(ns)*betas(ib)*betatypes(ibstsim,1))/2;
                    delta_lr_neg = 1 - (negpos_diff_fb(ilrf) + ...
                        phq(ns)*betas(ib)*betatypes(ibstsim,1))/2;
                    delta_lr_pos_conf = 1 + (negpos_diff_conf(ilrc) + ...
                        phq(ns)*betas(ib)*betatypes(ibstsim,2))/2;
                    delta_lr_neg_conf = 1 - (negpos_diff_conf(ilrc) + ...
                        phq(ns)*betas(ib)*betatypes(ibstsim,2))/2;

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
                        spe(ns,nb)= spe0(ns,nb+1) + ...
                            betas(ib)/10*betatypes(ibstsim,3); % with no learning to the next block
                    end

                    perf = mean(correct,3);

                end
                include_subj = all(perf>.6 & perf<=.85,2) & all(spe<1 & spe>0,2);
                % include_subj = true(1, ns);
                % sum(include_subj)

                % model_iter = model_iter+1;
                % fprintf('model iteration: %g\n', model_iter)

                spe = spe(include_subj,:);
                correct = correct(include_subj,:,:);
                conf = conf(include_subj,:,:);
                feedback = feedback(include_subj,:,:);
                task = groupTask(subj_group(include_subj),:);
                phq = phq(include_subj);
                fbblock = groupFeedback(subj_group(include_subj),:);
                group = subj_group(include_subj);

                for ibstrec = 1 :nbetatypes

                    fprintf('fb: %g,  conf: %g, beta: %g, bsim: %g, brec: %g\n', ...
                        ilrf, ilrc, ib, ibstsim, ibstrec)
                    fit = fit_globalSPE([], spe, correct, conf, feedback, ...
                        task, repmat(40, 1, 6), num2str(model_strs(ibstrec)), ...
                        njagiters, phq, tmpfolder);
                    dic_vals(ilrf,ilrc,ib,ibstsim,ibstrec) = fit.dic;

                    fprintf('fb: %g,  conf: %g, beta: %g, bsim: %g, brec: %g, dic: %g\n', ...
                        ilrf, ilrc, ib, ibstsim, ibstrec, fit.dic)

                    save('~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/data/model_selection_recovery3_3dif')

                end
            end
            initvals(4) = 1;
        end
        initvals(3) = 1;
    end
    initvals(2) = 1;
% end
