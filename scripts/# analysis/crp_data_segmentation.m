% for thresh_1 data
% script to segment crp data

path = getPath();
load([path, '\Neural data\stim_features\crp_feature_mat_thresh1.mat']);

% filter out 14 and 15/3
vault_exper = cellfun(@(faulty) cellfun(@isempty, strfind(table2save.Properties.RowNames,faulty)),{'14-Mar','15-Mar','13-Mar','16-Mar'},...
    'UniformOutput',false);
valid_ind = prod(cell2mat(vault_exper),2);
isnanCell = @(c) cell2mat(cellfun(@(x) sum(isnan(x)),c,'UniformOutput',false));
tbl = table2save(logical(valid_ind),:);
%% Change std to var

std_names = ~cellfun(@isempty, cellfun(@(c) strfind(c,'_std'),tbl.Properties.VariableNames, 'UniformOutput',false));
tbl(:,std_names) = num2cell(table2array(tbl(:,std_names)).^2);
%% Only mean mu?

tbl(:,~std_names) = [];
%% visualize using PCA

data = table2array(tbl);
names = cellfun(@(a) strjoin({a{2:end}},' '), cellfun(@(c) textscan(c,'%s', 'delimiter', '_'),tbl.Properties.VariableNames), 'UniformOutput',false);
figure()
boxplot(data,'orientation','horizontal','labels',names)

C = corr(data,data);

[pc,score,latent,tsquare,explained] = pca(table2array(tbl));
cumsum(explained)
figure()
h_pareto = pareto(explained);
set(h_pareto(2),'Color','k');
set(h_pareto(1),'FaceColor',[0.7294    0.8314    0.9569]);
xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('CRP data')
set(gca,'FontSize',14);
%%

figure()
scatter(score(:,1),score(:,2),'filled')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3nd Principal Component')
%%
pc2keep = {'amp std'};
a = cellfun(@(c) ismember(names,c),pc2keep, 'UniformOutput',false);
a = cellfun(@(c) strfind(names,c),pc2keep, 'UniformOutput',false);
pc_index = a{:};
if ~islogical( pc_index(1) )
    pc_index = ~cellfun(@(b) isempty(b), pc_index);
end
%%

num_pc = 2;
pc_index = (1:length(pc));
figure()
biplot(pc(:,1:num_pc),'scores',score(:,1:num_pc),'varlabels',{});
axis([-.26 0.6 -.51 .51]);
if num_pc >2 ; view([30 40]);end;
%% filter by...

trs_pc1 = 0.24;
trs1_idx = abs(pc(:,1)) >= trs_pc1;
figure(1);
ax1 = histogram(abs(pc(:,1)),16);
X1 = linspace(ax1.BinLimits(1),ax1.BinLimits(2),1e3);
Y1 = max(ax1.Values) .* (X1 >= trs_pc1);
hold on; 
ax12 = fill( [X1 fliplr(X1)],  [Y1 zeros(size(Y1))], 'r');
alpha(.2);
set(ax12,'EdgeColor',[0 0 0]);
set(ax12,'EdgeAlpha',0.5);
hold off;
title('Principal component 1 histogram (numBins = 16)')
xlabel('feature''s weight')
set(gca,'FontSize',20);
axis tight;

trs_pc2 = 0.2;
trs2_idx = abs(pc(:,2)) >= trs_pc2;
figure(2);histogram(abs(pc(:,2)),20)
ax2 = histogram(abs(pc(:,2)),16);
X1 = linspace(ax2.BinLimits(1),ax2.BinLimits(2),1e3);
Y1 = max(ax2.Values) .* (X1 >= trs_pc2);
hold on; 
ax12 = fill( [X1 fliplr(X1)],  [Y1 zeros(size(Y1))], 'r');
alpha(.2);
set(ax12,'EdgeColor',[0 0 0 ]);
set(ax12,'EdgeAlpha',0.5);
hold off;
axis tight;
title('Principal component 2 histogram (numBins = 16)')
xlabel('feature''s weight')
set(gca,'FontSize',20);

trs_r = 0.1;
trs_r_idx = sqrt(pc(:,2).^2 + pc(:,1).^2)  >= trs_r;
%% trs_1

num_pc = 2;
figure(11)
h_bi1 = biplot(pc(:,1:num_pc),'scores',score(:,1:num_pc),'varlabels',{});
line_list= findobj(h_bi1,'Tag','varline');
for k = 1:numel(line_list)
    set(line_list(k),'LineWidth',1.5,'Color',[0.8314    0.8157    0.7843]);
end
line_marker_list = findobj(h_bi1,'Tag','varmarker');
for k = 1:numel(line_marker_list)
    set(line_marker_list(k),'Marker','none','MarkerFaceColor',[ 0    0.4471    0.7412],'Color',[ 0    0.4471    0.7412]);
end
hold on;
h_bi2 = biplot(pc(trs1_idx,1:num_pc),'scores',score(:,1:num_pc),'varlabels',names(trs1_idx));
prettyBiplot( h_bi2 )
axis([-.26 0.6 -.51 .51]);
if num_pc >2 ; view([30 40]);end;
hold off;
font_size = 20;
xlabel(gca,'PC 1','FontSize',font_size)
ylabel(gca,'PC 2','FontSize',font_size)
title(gca,'PCA','FontSize',font_size);
%% trs_2

figure(12)
num_pc = 2;
h_bi1 = biplot(pc(:,1:num_pc),'scores',score(:,1:num_pc),'varlabels',{});
line_list= findobj(h_bi1,'Tag','varline');
for k = 1:numel(line_list)
    set(line_list(k),'LineWidth',1.5,'Color',[0.8314    0.8157    0.7843]);
end
line_marker_list = findobj(h_bi1,'Tag','varmarker');
for k = 1:numel(line_marker_list)
    set(line_marker_list(k),'Marker','none','MarkerFaceColor',[ 0    0.4471    0.7412],'Color',[ 0    0.4471    0.7412]);
end
hold on;
h_bi2 = biplot(pc(trs2_idx,1:num_pc),'scores',score(:,1:num_pc),'varlabels',names(trs2_idx));
prettyBiplot( h_bi2 )
axis([-.26 0.6 -.51 .51]);
if num_pc >2 ; view([30 40]);end;
hold off;
font_size = 20;
xlabel(gca,'PC 1','FontSize',font_size)
ylabel(gca,'PC 2','FontSize',font_size)
title(gca,'PCA','FontSize',font_size);
%%

names_temp = names;
names_temp(~pc_index) = {' '};
figure()
biplot(pc(:,1:3),'scores',score(:,1:3),'varlabels',names_temp(:));
axis([-.26 0.6 -.51 .51]);
view([30 40]);
%% plot explanation

% All nine variables are represented in this bi-plot by a vector, and the direction and length of the 
% vector indicate how each variable contributes to the two principal components in the plot. For example,
% the first principal component, on the horizontal axis, has positive coefficients for all nine variables. That
% is why the nine vectors are directed into the right half of the plot. The largest coefficients in the first principal
% component are the third and seventh elements, corresponding to the variables health and arts.

% The second principal component, on the vertical axis, has positive coefficients for the variables education, health,
% arts, and transportation, and negative coefficients for the remaining five variables. This indicates that the second component
% distinguishes among cities that have high values for the first set of variables and low for the second, and cities that have 
% the opposite.

% The variable labels in this figure are somewhat crowded. You can either exclude the VarLabels parameter when making the plot,
% or select and drag some of the labels to better positions using the Edit Plot tool from the figure window toolbar.

% This 2-D bi-plot also includes a point for each of the 329 observations, with coordinates indicating the score of each
% observation for the two principal components in the plot. For example, points near the left edge of this plot have the lowest
% scores for the first principal component. The points are scaled with respect to the maximum score value and maximum coefficient
% length, so only their relative locations can be determined from the plot.

% You can identify items in the plot by selecting Tools>Data Cursor from the figure window. By clicking a variable (vector), 
% you can read that variable's coefficients for each principal component. By clicking an observation (point), you can read that
% observation's scores for each principal component.