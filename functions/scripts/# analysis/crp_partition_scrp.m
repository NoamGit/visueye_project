% crp script for deciding of how many partitions should be used for the
% chirp and amplitude stimuli
load('C:\Users\noambox\Documents\Sync\Neural data\table_data\thres1\crp_table.mat');
load('C:\Users\noambox\Documents\Sync\Neural data\table_data\thres1\maps.mat');

% filter out 14 and 15/3
vault_exper = cellfun(@(faulty) cellfun(@isempty, strfind(crp_table.cell_id,faulty)),{'14-Mar','15-Mar','13-Mar','16-Mar'},...
    'UniformOutput',false);
valid_ind = prod(cell2mat(vault_exper),2);
isnanCell = @(c) cell2mat(cellfun(@(x) sum(isnan(x)),c,'UniformOutput',false));
crp_table = crp_table(logical(valid_ind),:);

data = crp_table.data;

% crop only chirp and amplitude stimulus times according to each cell
% stimuli
N = size(data,1);
count = 1;
[chirp_mat_S,chirp_mat_C, amp_mat_S, amp_mat_C] = deal({});
for k = 1:N
   key = cell2mat(crp_table(k,:).fkey);
   stim = stim_map(key).stim;
   i_crp = find(stim > 1.5 & stim < 2.5);
   i_amp = find(stim > 3 & stim < 4);
   for n = 1:size(stim_map(key).partition,2)
       i_part1 = find(i_crp >= stim_map(key).partition(1,n) & i_crp <= stim_map(key).partition(2,n));
       chirp_mat_S{count} =  data(k).S(i_crp(i_part1));
       chirp_mat_C{count} =  data(k).C(i_crp(i_part1));
       
       i_part2 = find(i_amp >= stim_map(key).partition(1,n) & i_amp <= stim_map(key).partition(2,n));
       amp_mat_S{count} = data(k).S(i_amp(i_part2));
       amp_mat_C{count} = data(k).C(i_amp(i_part2));
       count = count+1;
   end
end

%% for S data

% find maximum size and interp everything to this size
max_df_size = max(cellfun(@length, chirp_mat_S));
xq = linspace(0,max_df_size-1,max_df_size);
eqSize_crp = @(s) z1(reshape(interp1(linspace(0,length(s)-1,length(s)), s, xq,'spline'),[],1))';
chirp_mat_intrp = cellfun(eqSize_crp, chirp_mat_S, 'UniformOutput',false)';
figure(2); imagesc(cell2mat(chirp_mat_intrp))
title('Spike estimation to chirp section');xlabel('time lags');ylabel('trial');set(gca,'FontSize',11);

max_df_size = max(cellfun(@length, amp_mat_S));
xq = linspace(0,max_df_size-1,max_df_size);
eqSize_amp = @(s) z1(reshape(interp1(linspace(0,length(s)-1,length(s)), s, xq,'spline'),[],1))';
amp_mat_intrp = cellfun(eqSize_amp, amp_mat_S, 'UniformOutput',false)';
figure(3); imagesc(cell2mat(amp_mat_intrp))
title('Spike estimation to amplitude section');xlabel('time lags');ylabel('trial');set(gca,'FontSize',11);
%% for C data

max_df_size = max(cellfun(@length, chirp_mat_C));
xq = linspace(0,max_df_size-1,max_df_size);
chirp_mat_intrp = cellfun(eqSize_crp, chirp_mat_C, 'UniformOutput',false)';
figure(4); imagesc(cell2mat(chirp_mat_intrp))
title('Ca estimation to chirp section');xlabel('time lags');ylabel('trial');set(gca,'FontSize',11);

max_df_size = max(cellfun(@length, amp_mat_C));
xq = linspace(0,max_df_size-1,max_df_size);
amp_mat_intrp = cellfun(eqSize_amp, amp_mat_C, 'UniformOutput',false)';
figure(5); imagesc(cell2mat(amp_mat_intrp))
title('Ca estimation to amplitude section');xlabel('time lags');ylabel('trial');set(gca,'FontSize',11);
%% step detection - simple approach
chirp_mat_intrp = cellfun(eqSize_crp, chirp_mat_S, 'UniformOutput',false)';
amp_mat_intrp = cellfun(eqSize_amp, amp_mat_S, 'UniformOutput',false)';

M = count-1;
A_crp = cell2mat(chirp_mat_intrp);
bin_size = 3; % bin_size * 0.1 sec
filt = ones(M,bin_size);
B = conv2(A_crp, filt,'same');
mean_sample = M/2;
figure(6); imagesc(B);
mov_av_crp = B(round(mean_sample+1),:)./(M*bin_size);
figure(7);findpeaks(mov_av_crp,'MinPeakProminence',0.002);
title('Moving average on the chirp section');xlabel('time lags');ylabel(['average spike estimation binSize = ',num2str(bin_size)]);set(gca,'FontSize',11);
axis tight

A_amp = cell2mat(amp_mat_intrp);
bin_size = 6; % bin_size * 0.1 sec
filt = ones(M,bin_size);
B = conv2(A_amp, filt,'same');
mean_sample = M/2;
figure(8); imagesc(B);
mov_av_amp = B(round(mean_sample+1),:)./(M*bin_size);
figure(9);findpeaks(mov_av_amp,'MinPeakProminence',0.004);
%  plot(mov_av);
title('Moving average on the amp section');xlabel('time lags');ylabel(['average spike estimation binSize = ',num2str(bin_size)]);set(gca,'FontSize',11);
axis tight
%% plot with partition

amp_mat_intrp = cellfun(eqSize_amp, amp_mat_C, 'UniformOutput',false)';
chirp_mat_intrp = cellfun(eqSize_crp, chirp_mat_C, 'UniformOutput',false)';

A_amp = cell2mat(amp_mat_intrp);
A_crp = cell2mat(chirp_mat_intrp);

[~,crp_peaks] = findpeaks(mov_av_crp,'MinPeakProminence',0.002);
A_crp(:,crp_peaks(2:end-1)) = 0;
figure(10); imagesc(A_crp)
title('Spike estimation to chirp section');xlabel('time lags');ylabel('trial');set(gca,'FontSize',11);

[~,amp_peaks] = findpeaks(mov_av_amp,'MinPeakProminence',0.004);
A_amp(:,amp_peaks(2:end-1)) = 0;
figure(11); imagesc(A_amp)
title('Spike estimation to amp section');xlabel('time lags');ylabel('trial');set(gca,'FontSize',11);
