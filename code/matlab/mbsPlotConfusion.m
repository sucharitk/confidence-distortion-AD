
%%
load('data/model_selection_recovery_simfb.mat'); % with parameters -.04:.02:.04
dv1 = squeeze(dic_vals(:,:,:,1,:));
[mn,mnind1] = min(dv1,[],4);
mnind1(isnan(mn)) = NaN;

%%
load('data/model_selection_recovery_simconf.mat'); % with parameters -.04:.02:.04
dv2 = squeeze(dic_vals(:,:,:,2,:));
[mn,mnind2] = min(dv2,[],4);
mnind2(isnan(mn)) = NaN;

%%
load('data/model_selection_recovery_simaddbias.mat'); % with parameters -.04:.02:.04
dv3 = squeeze(dic_vals(:,:,:,3,:));
[mn,mnind3] = min(dv3,[],4);
mnind3(isnan(mn)) = NaN;

%% confusion matrix: p(fit model | simulated model)

cm = NaN(3);
sz = sum(~isnan(mnind1(:)));

cm(1,1) = sum(mnind1(:)==1)/sz;
cm(1,2) = sum(mnind1(:)==2)/sz;
cm(1,3) = sum(mnind1(:)==3)/sz;
cm(2,1) = sum(mnind2(:)==1)/sz;
cm(2,2) = sum(mnind2(:)==2)/sz;
cm(2,3) = sum(mnind2(:)==3)/sz;
cm(3,1) = sum(mnind3(:)==1)/sz;
cm(3,2) = sum(mnind3(:)==2)/sz;
cm(3,3) = sum(mnind3(:)==3)/sz;
cm

%% inversion matrix: p(simulated model | fit model)

mnind(1,:,:,:) = mnind1;
mnind(2,:,:,:) = mnind2;
mnind(3,:,:,:) = mnind3;

im = NaN(3);

im(1,1) = sum(mnind1(:)==1)/sum(mnind(:)==1);
im(1,2) = sum(mnind2(:)==1)/sum(mnind(:)==1);
im(1,3) = sum(mnind3(:)==1)/sum(mnind(:)==1);
im(2,1) = sum(mnind1(:)==2)/sum(mnind(:)==2);
im(2,2) = sum(mnind2(:)==2)/sum(mnind(:)==2);
im(2,3) = sum(mnind3(:)==2)/sum(mnind(:)==2);
im(3,1) = sum(mnind1(:)==3)/sum(mnind(:)==3);
im(3,2) = sum(mnind2(:)==3)/sum(mnind(:)==3);
im(3,3) = sum(mnind3(:)==3)/sum(mnind(:)==3);
im

%% 
rownames = {'Feedback dist.'; 'Confidence dist.'; 'Response bias'};
% imt = array2table(im, 'VariableNames',rownames);
% imt.Names = rownames;
figure
ax=heatmap(cm, 'Colormap', summer, 'FontSize', 15);
properties(ax)
ax.XDisplayLabels = rownames;
ax.YDisplayLabels = rownames;
ax.XLabel         = 'Simulated model';
ax.YLabel         = 'Fit model';
% ax.Title          = 'P(fit model | simulated model)';
colorbar;
clim([0 1])

%% 
rownames = {'Feedback dist.'; 'Confidence dist.'; 'Response bias'};
% imt = array2table(im, 'VariableNames',rownames);
% imt.Names = rownames;
figure
ax=heatmap(round(im*100)/100, 'Colormap', summer, 'FontSize', 15);
properties(ax)
ax.XDisplayLabels = rownames;
ax.YDisplayLabels = rownames;
ax.XLabel         = 'Simulated model';
ax.YLabel         = 'Fit model';
% ax.Title          = 'P(simulated model | fit model)';
colorbar;
clim([0 1])