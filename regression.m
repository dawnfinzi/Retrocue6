%%% computes regression for ERSP and WM accuracy (response error), computes 
%%% significance, and plots significance on the group level (using standardized 
%%% betaweights of the regression model) for the novel and standard conditions

% - load single trial ERSP matrices for each subject/IC
% - run regression (calls regressERSP_Z.m)
% - compute significance of standardardized beta coefficients over subjects
% - multiple comparison correction with fdr
% - mask standardized beta coefficients
% - plot average standardized beta coefficients for novel and standard
% - localize significant area of interest in novel trials and plot correlation 
% there between time x frequency power in that peak region and WM error
% (for visualization purposes only!)

close all

[folder, name, ext] = fileparts(which('regression.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

mypath = [folder, '/pp/STUDY/']; 

load([mypath 'Cluster_RF.mat'])
CLName = 'RF';
Cls_SUBJECT = RFCSJ;
Cls = RFCcomps;

for s = 1:length(Cls); 
    
    clear standard novel EEG_WM EEG_S EEG_N regNOVEL regSTAND
    
    %% reload individual datasets
    Participant_code = Cls_SUBJECT(s);
    str = char(Participant_code);
    subj_num = str(2:end); 
    if subj_num(1) == '0'
        subj_num = subj_num(2);
    end
    loadpath = [folder, '/pp/s', subj_num, '/'];

    EEG_WM = pop_loadset( [Cls_SUBJECT{s} '_WM_brain_ica.set'], loadpath);

    %% epoch datasets relative to the tone for retrocue
    EEG_S = pop_epoch( EEG_WM, {  '65530'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_S = eeg_checkset(EEG_S);

    EEG_N = pop_epoch( EEG_WM, {  '65529'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_N = eeg_checkset(EEG_N);

    %% create error vectors for each
    t = 1; 
    for i = 1:length(EEG_S.event)
        if EEG_S.event(i).type == '65530'
            standard(t) = EEG_S.event(i).error;
            t = t+1;
        end
    end

    t = 1; 
    for i = 1:length(EEG_N.event)
        if EEG_N.event(i).type == '65529'
            novel(t) = EEG_N.event(i).error;
            t = t+1;
        end
    end    

    %% load ERSP single trial file
    load([mypath 'ERSP_', CLName, '_TF_st_', num2str(s), '.mat']);
    
    %% run regression with time x frequency as independent var and WM error as dependent var
    regNOVEL = regressERSP_Z(tfdataNovel,novel, [-50 850], times, freqs);
    regSTAND = regressERSP_Z(tfdataStand,standard, [-50 850], times, freqs);
    
    % standardized beta coefficients
    b_vals_novel(s,:,:) = regNOVEL.b;
    b_vals_stand(s,:,:) = regSTAND.b;
end

save([mypath 'std_b_vals_', CLName, '.mat'], 'b_vals_novel', 'b_vals_stand'); 

%% compute significance against 0
zero_map = zeros(12,225,270);
zero_map = permute(zero_map, [2 3 1]); 

b_v_n= permute(b_vals_novel, [2 3 1]);
b_v_s= permute(b_vals_stand, [2 3 1]);

[stats, df, pNovel, surrog] = statcond( {b_v_n, zero_map}, 'mode', 'bootstrap', 'naccu', 1000);
[stats, df, pStand, surrog] = statcond( {b_v_s, zero_map}, 'mode', 'bootstrap', 'naccu', 1000);

%% PLOTTING
figpath = [folder, '/figures/']; 

%% which timewindow to plot
lat = [-50 850];
for i = 1:length(lat)
[minDifferenceValue, indexLat(i)] = min(abs(times - lat(i)));
end

%% up to which frequency plot the ERSP image
freqind = find(freqs == 30);
freqind1 = find(freqs == 4);

%% correct for multiple comparisons
[pStand1,pmaskedStand] = fdr(pStand(freqind1:freqind,:),0.05, 'Parametric');
[pNovel1,pmaskedNovel] = fdr(pNovel(freqind1:freqind,:),0.05, 'Parametric');

%% average beta coefficients for plotting and limit frequency range
b_vals_n = squeeze(mean(b_vals_novel));
b_vals_s = squeeze(mean(b_vals_stand));

b_vals_n = b_vals_n(freqind1:freqind,:);
b_vals_s = b_vals_s(freqind1:freqind,:);

%% latencies at which to plot vertical lines
Novel_newLat = [0 800]; % when target appears

clear indexNOVEL
clear indexGo
clear indexSTOP
clear indexFAIL

for i = 1:length(Novel_newLat)
[minDifferenceValue, indexNOVEL(i)] = min(abs(times - Novel_newLat(i)));
end

%% reshape data for plotting
AVnovel=reshape(b_vals_n,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskednovel = reshape(pmaskedNovel, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

AVstand=reshape(b_vals_s,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskedstand = reshape(pmaskedStand, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

%% set maximum for scaling image
MAXN = .5; 

fu1 = figure ( ); 
fu2 = figure ( );

%% plot novel - avg std. beta coefficients masked by significance
figure ( fu1 );
stat.dimord = 'chan_freq_time';
stat.freq = freqs(freqind1:freqind);
stat.time = times(indexLat(1):indexLat(2));
stat.cfg.channel = {'Chosen'};
stat.label = {'Chosen'}; 
stat.mask = maskednovel;
stat.effect= AVnovel;
stat.stat = AVnovel;
stat.cfg.equency=[];
j=jet;
colormap(j)
cfg = [];
cfg.channel       = {'Chosen'};
cfg.ylim= [freqs(freqind1) freqs(freqind)];
cfg.zlim           = [-MAXN MAXN];
cfg.maskparameter = 'mask';       % use significance to mask the power
cfg.maskalpha     = 0.3;   
cfg.renderer      = 'openGL';     % painters does not support opacity, openGL does
cfg.colorbar      = 'yes';
cfg.parameter     = 'effect';     % display the power        % make non-significant regions 30% visible
ft_singleplotTFR(cfg, stat);
title('significant power changes (p<0.05, corrected)')
colorbar
cbar_axes = colorbar('peer',gca);
get(cbar_axes)
set(get(cbar_axes,'ylabel'),'String', 'dB', 'fontsize',16);
set(gca,'fontsize',30)
hold on
plot(get(gca,'xlim'), [8 8], 'Color','k', 'LineWidth',1, 'LineStyle', '--'); % Adapts to x limits of current axes
hold on
plot(get(gca,'xlim'), [12 12], 'Color','k', 'LineWidth',1, 'LineStyle', '--'); % Adapts to x limits of current axes
hold on
plot(get(gca,'xlim'), [20 20], 'Color','k', 'LineWidth',1, 'LineStyle', '--'); % Adapts to x limits of current axes
hold on
for k = 1:length(Novel_newLat)
line([[times([indexNOVEL(k)])] [times([indexNOVEL(k)])]],[freqs(freqind1) freqs(freqind)], 'LineWidth',2);
end
text([[times([indexNOVEL(1)])+15] [times([indexNOVEL(1)])+15]], [28 28],'Novel', 'fontsize',30);
text([[times([indexNOVEL(2)])+15] [times([indexNOVEL(2)])+15]], [26 26],'Target', 'fontsize',30);
ylabel('frequency (Hz)','fontsize',30)
xlabel('time (ms)','fontsize',30)
title('Novel tone','fontsize',30)
set(gca,'fontsize',30)

saveas(fu1, [figpath 'WM_regZ_AV_Novel_', CLName, '.fig'])
print(fu1,  '-depsc2 ', [figpath 'WM_regZ_AV_Novel_', CLName])


%% plot standard - avg std. beta coefficients masked by significance
figure ( fu2 );
stat.dimord = 'chan_freq_time';
stat.freq = freqs(freqind1:freqind);
stat.time = times(indexLat(1):indexLat(2));
stat.cfg.channel = {'Chosen'};
stat.label = {'Chosen'}; 
stat.mask = maskedstand;
stat.effect= AVstand;
stat.stat = AVstand;
stat.cfg.equency=[];
j=jet;
colormap(j)
cfg = [];
cfg.channel       = {'Chosen'};
cfg.ylim= [freqs(freqind1) freqs(freqind)];
cfg.zlim           = [-MAXN MAXN];
cfg.maskparameter = 'mask';       % use significance to mask the power
cfg.maskalpha     = 0.3;   
cfg.renderer      = 'openGL';     % painters does not support opacity, openGL does
cfg.colorbar      = 'yes';
cfg.parameter     = 'effect';     % display the power        % make non-significant regions 30% visible
ft_singleplotTFR(cfg, stat);
title('significant power changes (p<0.05, corrected)')
colorbar
cbar_axes = colorbar('peer',gca);
get(cbar_axes)
set(get(cbar_axes,'ylabel'),'String', 'dB', 'fontsize',16);
set(gca,'fontsize',30)
hold on
plot(get(gca,'xlim'), [8 8], 'Color','k', 'LineWidth',1, 'LineStyle', '--'); % Adapts to x limits of current axes
hold on
plot(get(gca,'xlim'), [12 12], 'Color','k', 'LineWidth',1, 'LineStyle', '--'); % Adapts to x limits of current axes
hold on
plot(get(gca,'xlim'), [20 20], 'Color','k', 'LineWidth',1, 'LineStyle', '--'); % Adapts to x limits of current axes
hold on
for k = 1:length(Novel_newLat)
line([[times([indexNOVEL(k)])] [times([indexNOVEL(k)])]],[freqs(freqind1) freqs(freqind)], 'LineWidth',2);
end
text([[times([indexNOVEL(1)])+15] [times([indexNOVEL(1)])+15]], [28 28],'Standard', 'fontsize',30);
text([[times([indexNOVEL(2)])+15] [times([indexNOVEL(2)])+15]], [26 26],'Target', 'fontsize',30);
ylabel('frequency (Hz)','fontsize',30)
xlabel('time (ms)','fontsize',30)
title('Standard tone','fontsize',30)
set(gca,'fontsize',30)

saveas(fu2, [figpath 'WM_regZ_AV_Stand_', CLName, '.fig'])
print(fu2,  '-depsc2 ', [figpath 'WM_regZ_AV_Stand_',CLName])

%% FOR VISUALIZATION PURPOSES ONLY

%% Single trial - all subjects plot (power from significant cluster)
 
for s = 1:length(Cls); 
    
    clear standard novel EEG_WM EEG_S EEG_N tfmaptmp tfmapsel testing power_comp power_avg minMask*
    
    %% reload individual datasets
    Participant_code = Cls_SUBJECT(s);
    str = char(Participant_code);
    subj_num = str(2:end); 
    if subj_num(1) == '0'
        subj_num = subj_num(2);
    end
    loadpath = [folder, '/pp/s', subj_num, '/'];
    
    EEG_WM = pop_loadset( [Cls_SUBJECT{s} '_WM_brain_ica.set'], loadpath);

    %% epoch datasets relative to the tone for retrocue
    EEG_S = pop_epoch( EEG_WM, {  '65530'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_S = eeg_checkset(EEG_S);

    EEG_N = pop_epoch( EEG_WM, {  '65529'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_N = eeg_checkset(EEG_N);

    %% create error vectors for each
    t = 1; 
    for i = 1:length(EEG_S.event)
        if EEG_S.event(i).type == '65530'
            standard(t) = EEG_S.event(i).error;
            t = t+1;
        end
    end

    t = 1; 
    for i = 1:length(EEG_N.event)
        if EEG_N.event(i).type == '65529'
            novel(t) = EEG_N.event(i).error;
            t = t+1;
        end
    end    

    %% load ERSP single trial file
    load([mypath 'ERSP_', CLName, '_TF_st_', num2str(s), '.mat']);
    
    %% find the power (normalized) at only the relevant peak time window
    timewindow = [-50 850];
    tfmaptmp = abs(tfdataNovel); 
    baseind = find(times < 0);
    timesind = find(times > timewindow(1) & times < timewindow(2));
    tfmapsel = tfmaptmp(:,timesind,:);
    
    baseline = repmat(mean(tfmaptmp(:,baseind,:),2),1,size(tfmapsel,2),1);
    Zmap = log10(tfmapsel)./log10(baseline);
    Zmap = Zmap(freqind1:freqind,:,:);
    
    [minValue,minIdx]=min(pNovel(:));
    minMask = pNovel <= minValue;
    minMask1 = repmat(minMask(freqind1:freqind,:),[1,1,size(tfmaptmp,3)]);
    testing = Zmap.*minMask1; 

    freqind2 = find(freqs == 20);
    freqind3 = find(freqs == 15);
    power_comp = testing(freqind3:freqind2,:,:); % take only freq range 15 to 20 to localise power 

    for i = 1:size(power_comp,3)
        clear count total
        count = sum(sum(power_comp(:,:,i)>0));
        total = sum(sum(power_comp(:,:,i)));

        power_avg(i) = total/count;
    end

    compiled.response_error{s} = zscore(novel); % response error normalized for each trial
    compiled.re_basic{s} = novel; % response error for each trial
    compiled.power{s} = power_avg; % average power in peak time x frequency window for each trial

end
%% plot correlations divided by subject 

fu3 = figure ( );
figure ( fu3 );
for i = 1:12
    subplot(3,4,i); 
    x = compiled.response_error{i}/46.6;
    y = compiled.power{i};
    scatter(x,y, 'r','LineWidth', 1);
    set(gca,'xlim', [min(x)-.01 max(x)+.01]);

    hold  on
    x_lim = [min(x)-.01 max(x)+.01];
    coef_fit = polyfit(x,y,1);
    y_fit = polyval(coef_fit,x_lim);
    plot(x_lim,y_fit,'k', 'LineWidth',1.5);
end

saveas(fu3, [figpath 'reg_corr_vis_', CLName, '.fig'])
print(fu3,  '-depsc2 ', [figpath 'reg_corr_vis_', CLName])

%% plot correlations for each trial and each subject on one plot
clear x y
x = []; y = [];
fu4 = figure ( );
figure ( fu4 );

for i = 1:12
    x = [x compiled.response_error{i}/46.6];
    y = [y compiled.power{i}];
end

scatter(x,y,[],[.1 .1 .1],'LineWidth', 1);
set(gca,'xlim', [min(x)-.01 max(x)+.01]);

hold  on
x_lim = [min(x)-.01 max(x)+.01];
coef_fit = polyfit(x,y,1);
y_fit = polyval(coef_fit,x_lim);
plot(x_lim,y_fit,'k', 'LineWidth',1.5);

saveas(fu4, [figpath 'reg_corr_overall_', CLName, '.fig'])
print(fu4,  '-depsc2 ', [figpath 'reg_corr_overall_', CLName])