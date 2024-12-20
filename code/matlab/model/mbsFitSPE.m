function mbsData = mbsFitSPE(mbsData, modelName, nSamples, ...
    zconf, mhqNum, useReportedSPEPrior)
%
% wrapper function to call the function to fit the model
%
conf = mbsData.conf;

if zconf
    % z transform confidence data
    nanconf = isnan(conf);
    reconf = reshape(conf, [size(conf,1) size(conf,2)*size(conf,3)]);
    renanconf = reshape(nanconf, [size(conf,1) size(conf,2)*size(conf,3)]);
    nancols = all(renanconf);
    reconf(:,~nancols) = zscore(reconf(:,~nancols), [], 2);

    reconf = reconf - repmat(min(reconf, [], 2), [1 size(reconf,2)]); % min to 0
    reconf = reconf ./ repmat(max(reconf, [], 2), [1 size(reconf,2)]); % max to 1
    conf = reshape(reconf, [size(conf,1) size(conf,2) size(conf,3)]);
end

corr = mbsData.corr;
feedback = mbsData.feedback;
spe = mbsData.spe;
task = mbsData.groupTask;
fbBlock = mbsData.fbBlock;
fbBlock(fbBlock==2) = -1;

modelNumber = modelName(end-1:end);
if modelNumber(1)=='_', modelNumber=modelNumber(end); end
modelNum = str2num(modelNumber);

switch mbsData.expNum
    case 1
        ntrials = repmat(40, 1, 6);
        if modelNum<5
            % models 0-4 for comparison
            covariate = mbsData.questAvg(:,1);

            if zconf
                fitName = 'fitZ';
            else
                fitName = 'fit';
            end

        else
            % fit depression scores only with a diff in beta feedback model
            % (no conf asymmetry) - to compare with a posterior bias
            % model
            covariate = mbsData.questAvg(:,mhqNum);
            if zconf
                fitName = 'fitRegZ';
            else
                fitName = 'fitReg';
            end
          
        end


    case 2
        ntrials = repmat(40, 1, 6);
        ntrials([4 6]) = 20;

        if modelNum<5
            % models 0-4 for comparison
            covariate = mbsData.questAvg(:,1);
            if zconf
                fitName = 'fitZ';
            else
                fitName = 'fit';
            end

        elseif modelNum==11
            % pass all 3 transdiag axes to model 11
            covariate = mbsData.questAvg;
            if zconf
                fitName = 'fitRegZ';
            else
                fitName = 'fitReg';
            end

        else
            % fit depression scores only with a diff in beta feedback model
            % (no conf asymmetry) - to compare with a posterior bias
            % model
            covariate = mbsData.questAvg(:,mhqNum);
            if zconf
                fitName = 'fitRegZ';
            else
                fitName = 'fitReg';
            end

        end

end


if any(covariate<0)
    % for z-scored values make the covariate non-negative
    covariate = (covariate - min(covariate));
end


% subj 202 and 208 for run 4 gave conf=1 for all trials, so fit is crashing - make
% their conf values .999
remSubj = [];

switch mbsData.expNum
    case 1
        % conf(183,4,:) = .999;
        % conf(189,4,:) = .999;
        conf(202,4,:) = .999;
        conf(208,4,:) = .999;
        useReportedSPEPrior = false;
        % remSubj = [202,208];
    case 2
        conf(144,4,1:20) = .999;
        useReportedSPEPrior = false;
        % remSubj = [144];
end

spe(remSubj,:) = [];
conf(remSubj,:,:) = [];
corr(remSubj,:,:) = [];
feedback(remSubj,:,:) = [];
task(remSubj,:) = [];
covariate(remSubj) = [];
fbBlock(remSubj,:) = [];

% conf(conf==1) = .99;
% conf(conf==0) = .01;

if useReportedSPEPrior
    spe0 = mbsData.spe0;
else
    spe0 = [];
end

tmpfolder = 'tmpjags4';

fit = fit_globalSPE(spe0,...
    spe,corr,...
    conf,feedback,...
    task, ntrials,...
    modelNumber, nSamples, covariate, tmpfolder, fbBlock);

modelNameFull = [modelName '_q' num2str(mhqNum(end))];

mbsData.(fitName).(modelNameFull).spe_est = fit.spe_est;
mbsData.(fitName).(modelNameFull).dic = fit.dic;
if isfield(fit, 'samples')
    mbsData.(fitName).(modelNameFull).samples = fit.samples;
end
if isfield(fit, 'stats')
    mbsData.(fitName).(modelNameFull).stats = fit.stats;
end

end
