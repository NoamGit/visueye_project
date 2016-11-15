% feature analysis script
% data structure is [qi, mi, csd]

cd('C:\Users\noambox\Documents\Sync');
load('.\Neural data\stat_features\all_feat_spectral_last.mat');
load('.\Neural data\ori_table.mat');
% global ori_features_filt
global lbl
global crp_features_filt
global crp_table
global crp_features
%% annotate each data point with a label

fig = figure;
crp_features_filt = crp_features(crp_features.label == 1,:);
lbl = [crp_features_filt.datel ,crp_features_filt.filename];
l = length(crp_features_filt.QI);
c = (1+z1( crp_features_filt.QI )) .*255;
h0 = scatter3(crp_features_filt.QI,crp_features_filt.MI,crp_features_filt.CSD, 30,c);
h0.MarkerFaceColor = [0 .70 .70];

xlabel('QI');
ylabel('MI');
zlabel('CSD');
title('CRP features')
set(gca,'FontSize',11);
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@crpInteractive)