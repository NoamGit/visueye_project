function [ feature_struct ] = cnt_feature_extraction( cell_data, h )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% get summary features for cell
[ mean_response, keys_cpd, max_contrast,cell_data, group_indx ] = cntSummary(h, cell_data); 

% ANOVA testing for significane
group_cpd = group_indx(:,2);
[p,t,stats] = anova1(cell_data.data.S,group_cpd,'off');
[c,m,~,nms] = multcompare(stats,'Alpha',0.05);

end

