%% script - test spike inference MCMC

[fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
%% load data

sheet_cell = xlsread([pathName fName]);
time = sheet_cell(:,1);
dt = double(time(2)-time(1));
fs = 1/dt;
off_period = 15*fs+1;
stim = sheet_cell(:,2);
str = textscan(fName,'%s','delimiter','_');
full_field_s = [{'CRP'},{'WGN'}];
full_field = 0;
if(any(ismember(full_field_s,str{1})))
    sig_mat = sheet_cell(off_period:end,(3:end-1)); 
    dr = sheet_cell(off_period:end,end); % dark region in the FOV
    full_field = 1;
else
    sig_mat = sheet_cell(off_period:end,(3:end)); 
end
num_sig = size(sig_mat,2);
time = time(1:size(sig_mat,1));
%% load artifact

[fName, pathName] = uigetfile('D:\# Projects (Noam)\# SLITE\# DATA\*.xlsx','Select Cell''s Data Sheet');
sheet = xlsread([pathName fName]);
arti_mat = sheet(:,(3:end)); 
%% DF (data is after dark noise removal)
    
Fs_filt = 10;  % Sampling Frequency
N     = 5;     % Order
Fstop = 0.02;  % Stopband Frequency
Astop = 80;    % Stopband Attenuation (dB)
h = fdesign.lowpass('n,fst,ast', N, Fstop, Astop, Fs_filt);
Hd = design(h, 'cheby2', 'SystemObject', true);

for k = 1:num_sig
    
    x = sig_mat(:,k);
    
    % step 1 - removing low frq trend with LPF
    trend(:,k) = filtfilt(Hd.SOSMatrix,Hd.ScaleValues,x);
    x_detrend = x - trend(:,k); % x1 is the detrended signal
    
    % step 2 -  cut 15 first sec from data and artifact
    x_detrend = x_detrend(15*fs:end);
    x_detrend = x_detrend + mean(x); % avoiding deviding by 0 with bl
    t_trim = time(15*fs:end);
    
    % step 3 - Baseline estimation according lowest 15%
    Delta_t = 50 * fs; % sec * sr
    s_size = size(x_detrend,1);
    T = t_trim(end);
    num_bin = floor( s_size / Delta_t );
    s_mat_binned = reshape(x_detrend(1:num_bin*Delta_t,:),Delta_t, num_bin, []);
    t_mat_binned = reshape(t_trim(1:num_bin*Delta_t,:),Delta_t, num_bin, []);
    perc = .15; % the percent that the hist values under are probabely bsl
    quantile_vec = quantile(s_mat_binned, perc, 1);
    indx_mat = bsxfun(@(a,b) a < b,s_mat_binned,quantile_vec);
    Y_val = s_mat_binned(indx_mat); X_val = t_mat_binned(indx_mat);
    p = polyfit( X_val, Y_val,1);
    bl(:,k) = polyval(p,t_trim);
    df_mat(:,k) = (x_detrend-bl(:,k))./abs(bl(:,k));
end
%% plotting test

k = 4;
figure(6)
subplot(311);
plot([norm_nc(sig_mat(15*fs:end,k),5) norm_nc(trend(15*fs:end,k),5)]);axis tight
subplot(312);
plot([norm_nc(sig_mat(15*fs:end,k)-trend(15*fs:end,k),5) norm_nc(bl(:,k),5)]);axis tight
subplot(313);
plot(df_mat(:,k),'r');axis tight
%% Artifact removal
% function input order(signal, arti, on_time, filt_bs, safety_offset, filt_thres ,fs, maxl_phase)

[ sm_proc, a_proc ] = ArtiRemOrientation(df_mat, arti_mat(15*fs:end,:), 5, 11, 30, 1e3, 10, 16);
sm_proc_norm = arrayfun(@(x) norm_nc(sm_proc(:,x)), (1:num_sig),'UniformOutput',false);
figure(7); imagesc(cell2mat(sm_proc_norm)')
%%

k = 2;
figure(6)
subplot(311);
plot([norm_nc(sig_mat(15*fs:end,k),5) norm_nc(trend(15*fs:end,k),5)]);axis tight
subplot(312);
plot([norm_nc(sig_mat(15*fs:end,k)-trend(15*fs:end,k),5) norm_nc(bl(:,k),5)]);axis tight
subplot(313);
plot(sm_proc(:,k),'r');axis tight
%% MCMC

k = 4; % 34, 10, 2, 4, 9 ,48
x = sm_proc(:,k);

% initialize params
params.B = 0;            % number of burn in samples (default 200)
params.Nsamples = 100;    % number of samples after burn in (default 500)
params.f = 0.1;          % time step size
% params.g = exp(-params.f./9);
params.p = 2;                 % order of AR model (p == 1 or p == 2, default 1)
params.marg = 0;         % flag for marginalized sampler (default 0)
% params.upd_gam = 1;      % flag for updating gamma (default 0)

SAMPLES = cont_ca_sampler(z1(x),params);
C = make_mean_sample(SAMPLES,z1(x));
S = mean(samples_cell2mat(SAMPLES.ss,size(C,2)));
YrA = z1(x) - C';
P.b = mean(SAMPLES.Cb);
P.c1 = mean(SAMPLES.Cin);
P.neuron_sn = sqrt(mean(SAMPLES.sn2));
P.gn = mean(exp(-params.f./SAMPLES.g));
P.samples_mcmc = SAMPLES; % FN added, a useful parameter to have.
%% Plot 

Pl.n    = S; 
Pl.n(Pl.n==0)=NaN; % true spike train (0's are NaN's so they don't plot)
Pl.lw   = 1;                    % line width
figure(8);

subplot(211);
hold on;
h3 = plot(t_trim,(x),'k','LineWidth',Pl.lw);
h3.Color(4)=0.8;
h4 = plot(t_trim,C,'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
    'LineWidth',Pl.lw);
h4.Color(4)=0.5;
set(gca,'FontSize',14,'Box','off');
xlabel('Time [sec]');
ylabel('\Delta{F}/F');
axis('tight')
hold off;

subplot(212);
hold on; 
plot(t_trim,norm_nc(x,1),'Color',[0.709215695858002 0.909803926944733 0.850980401039124]);
h2 = plot(t_trim,norm_nc(C,5),'Color',[0.584313750267029 0.388235300779343 0.388235300779343],...
    'LineWidth',Pl.lw);
h2.Color(4)=0.5;
disp(['tau is - ',num2str(mean(-dt./log(P.gn))),' [sec]']);
disp(['SSE is - ',num2str(sum((norm_nc(x,1)' - norm_nc(C,5)).^2))]);
h1 = stem(t_trim,-1*Pl.n,'Color',[0.494117647409439 0.494117647409439 0.494117647409439]);
h1.ShowBaseLine ='off';
h1.Color(4)=0.01;
set(gca,'FontSize',14,'Box','off');
xlabel('Time [sec]');
ylabel('normelized C');
legend('F','est. C','est n_t');
set(legend,'Location','northwestoutside');
axis('tight')
hold off
%%   Variables:
%   y:      raw fluorescence data (vector of length(T))
%   c:      denoised calcium concentration (Tx1 vector)
%   b:      baseline concentration (scalar)
%  c1:      initial concentration (scalar)
%   g:      discrete time constant(s) (scalar or 2x1 vector)
%  sn:      noise standard deviation (scalar)
%  sp:      spike vector (Tx1 vector)
%% params.g              discrete time constant(s) (estimated if not provided)
% params.sn             initializer for noise (estimated if not provided)
% params.b              initializer for baseline (estimated if not provided)
% params.c1             initializer for initial concentration (estimated if not provided)
% params.marg           flag for marginalized sampler (default 0)
% params.upd_gam        flag for updating gamma (default 0)
% params.gam_step       number of samples after which gamma is updated (default 50)
% params.std_move       standard deviation of shifting kernel (default 3*Dt)
% params.add_move       number of add moves per iteration (default T/100)
% params.init           initial sample 
% params.f              imaging rate (default 1)
% params.p              order of AR model (p == 1 or p == 2, default 1)
% params.defg           default discrete time constants in case constrained_foopsi cannot find stable estimates
% params.TauStd         standard deviation for time constants in continuous time (default [0.2,2])
