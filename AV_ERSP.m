%%% average cluster ERSP images for each condition, compute significance,
%%% and plot differences for Novel-Standard and Successful-Failed Stop
%%% trials

% - load ERSP matrices
% - compute common baseline
% - compute significance relative to baseline over subjects
% - compute significance and in between conditions for Novel vs Standard and Successful stop vs fail
% - multiple comparison correction with fdr
% - mask ERSPs
% - plot average ERSPs
% - plot difference between Novel-Standard and Successful-Failed Stop
 
clear all
close all
clc

[folder, name, ext] = fileparts(which('AV_ERSP.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

mypath = [folder '/pp/STUDY/'];

CLName = 'RF';
load([mypath 'ERSP_RF_TF.mat']);
load([mypath 'Cluster_RF.mat']);

Cls_SUBJECT = RFCSJ;
Cls = RFCcomps;

mkdir(folder,'figures');
figpath = [folder '/figures/'];

% rename for consistency
ERSPStop = ERSPSucc;
BASEStop = BASESucc;

%% which timewindow to plot
lat = [-50 850];
for i = 1:length(lat)
[minDifferenceValue, indexLat(i)] = min(abs(times - lat(i)));
end

%% up to which frequency plot the ERSP image
freqind = find(freqs == 30);
freqind1 = find(freqs == 4);

%% baseline prep for computing significance

Bnovel = repmat(BASENovel(:,freqind1:freqind),[1,1, length(times(indexLat(1):indexLat(2)))]);
Bstand = repmat(BASEStand(:,freqind1:freqind),[1,1, length(times(indexLat(1):indexLat(2)))]);

Bfail = repmat(BASEFail(:,freqind1:freqind),[1,1, length(times(indexLat(1):indexLat(2)))]);
Bgo = repmat(BASEGo(:,freqind1:freqind),[1,1, length(times(indexLat(1):indexLat(2)))]);
Bstop = repmat(BASEStop(:,freqind1:freqind),[1,1, length(times(indexLat(1):indexLat(2)))]);

%% compute common baseline for novel and standards
Bnovelstand = (BASENovel(:,freqind1:freqind) + BASEStand(:,freqind1:freqind))/2;
Bnovelstand_diff = repmat(Bnovelstand,[1,1,length(times(indexLat(1):indexLat(2)))]);

%% compute common baseline for stop and fail trials
Bgostop = (BASEStop(:,freqind1:freqind) + BASEFail(:,freqind1:freqind))/2;
Bgostop_diff = repmat(Bgostop,[1,1,length(times(indexLat(1):indexLat(2)))]);

%% compile data for sig testing
Abs_novel = (ERSPNovel(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bnovel(:,:,:));
Abs_stand = (ERSPStand(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bstand(:,:,:));

Abs_go = (ERSPGo(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bgo(:,:,:));
Abs_fail = (ERSPFail(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bfail(:,:,:));
Abs_stop = (ERSPStop(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bstop(:,:,:));

%% reshape data
Abs_novel = permute(Abs_novel, [2 3 1]);
BNovel = permute(Bnovel, [2 3 1]);
Abs_stand = permute(Abs_stand, [2 3 1]);
BStand = permute(Bstand, [2 3 1]);

Abs_fail = permute(Abs_fail, [2 3 1]);
BFail = permute(Bfail, [2 3 1]);
Abs_go = permute(Abs_go, [2 3 1]);
BGo = permute(Bgo, [2 3 1]);
Abs_stop = permute(Abs_stop, [2 3 1]);
BStop = permute(Bstop, [2 3 1]);

Novel = {Abs_novel BNovel};
Stand = {Abs_stand BStand};

Fail = {Abs_fail BFail};
Go = {Abs_go BGo};
Stop = {Abs_stop BStop};

%% computing significance relative to baseline 
[stats, df, pStand, surrog] = statcond( Stand, 'mode', 'bootstrap', 'naccu', 2000);
[stats, df, pNovel, surrog] = statcond( Novel, 'mode', 'bootstrap', 'naccu', 2000);

[stats, df, pGo, surrog] = statcond( Go, 'mode', 'bootstrap', 'naccu', 2000);
[stats, df, pFail, surrog] = statcond( Fail, 'mode', 'bootstrap', 'naccu', 2000);
[stats, df, pStop, surrog] = statcond( Stop, 'mode', 'bootstrap', 'naccu', 2000);

%% correct with fdr
[pStand,pmaskedStand] = fdr(pStand,0.05, 'Parametric');
[pNovel,pmaskedNovel] = fdr(pNovel,0.05, 'Parametric');

[pGo,pmaskedGo] = fdr(pGo,0.05, 'Parametric');
[pFail,pmaskedFail] = fdr(pFail,0.05, 'Parametric');
[pStop,pmaskedStop] = fdr(pStop,0.05, 'Parametric');

%% get average values for plotting under significance mask
AV_Novel = squeeze(mean(ERSPNovel(:,freqind1:freqind,[indexLat(1):indexLat(2)]), 1));
AV_Stand = squeeze(mean(ERSPStand(:,freqind1:freqind,[indexLat(1):indexLat(2)]), 1));

AV_Go = squeeze(mean(ERSPGo(:,freqind1:freqind,[indexLat(1):indexLat(2)]), 1));
AV_Fail = squeeze(mean(ERSPFail(:,freqind1:freqind,[indexLat(1):indexLat(2)]), 1));
AV_Stop = squeeze(mean(ERSPStop(:,freqind1:freqind,[indexLat(1):indexLat(2)]), 1));

%% mask significance
maskNovel = (AV_Novel.*pmaskedNovel);
maskStand = (AV_Stand.*pmaskedStand);

maskGo = (AV_Go.*pmaskedGo);
maskFail = (AV_Fail.*pmaskedFail);
maskStop = (AV_Stop.*pmaskedStop);

%% determine latencies at which to plot vertical lines
Novel_newLat = [0 800]; % when target appears
GO_newLat = [0 536]; % average response on go trials
STOP_newLat = [0 252 540]; % average SSD and then avg SSRT
FAILs_newLat = [0 252 454]; % average SSD and then avg response on failed trials

clear indexNOVEL
clear indexGo
clear indexSTOP
clear indexFAIL

for i = 1:length(Novel_newLat)
[minDifferenceValue, indexNOVEL(i)] = min(abs(times - Novel_newLat(i)));
end

for i = 1:length(GO_newLat)
[minDifferenceValue, indexGo(i)] = min(abs(times - GO_newLat(i)));
end

for i = 1:length(STOP_newLat)
[minDifferenceValue, indexSTOP(i)] = min(abs(times - STOP_newLat(i)));
end

for i = 1:length(STOP_newLat)
[minDifferenceValue, indexFAIL(i)] = min(abs(times - FAILs_newLat(i)));
end

%% reshape data for plotting
AVnovel=reshape(AV_Novel,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskednovel = reshape(pmaskedNovel, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

AVstand=reshape(AV_Stand,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskedstand = reshape(pmaskedStand, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

AVstop=reshape(AV_Stop,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskedstop = reshape(pmaskedStop, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

AVfail=reshape(AV_Fail,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskedfail = reshape(pmaskedFail, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

AVgo=reshape(AV_Go,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskedgo = reshape(pmaskedGo, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

%% sey maximum for scaling image

MAXS = 1.5; MAXN = 1;

fu1 = figure ( ); 
fu2 = figure ( );
fu3 = figure ( ); 
fu4 = figure ( );
fu5 = figure ( );

%% plot Novel
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
 saveas(fu1, [figpath 'WM_' CLName '_AV_Novel.fig'])
 print(fu1,  '-depsc2 ', [figpath 'WM_' CLName '_AV_Novel'])


%% plot standard
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
 saveas(fu2, [figpath 'WM_' CLName '_AV_Stand.fig'])
 print(fu2,  '-depsc2 ', [figpath 'WM_' CLName '_AV_Stand'])


%% plot stop
figure ( fu3 );
stat.dimord = 'chan_freq_time';
stat.freq = freqs(freqind1:freqind);
stat.time = times(indexLat(1):indexLat(2));
stat.cfg.channel = {'Chosen'};
stat.label = {'Chosen'}; 
stat.mask = maskedstop;
stat.effect= AVstop;
stat.stat = AVstop;
stat.cfg.equency=[];
j=jet;
colormap(j)
cfg = [];
cfg.channel       = {'Chosen'};
cfg.ylim= [freqs(freqind1) freqs(freqind)];
cfg.zlim           = [-MAXS MAXS];
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
for k = 1:length(STOP_newLat)
line([[times([indexSTOP(k)])] [times([indexSTOP(k)])]],[freqs(freqind1) freqs(freqind)], 'LineWidth',2);
end
text([[times([indexSTOP(1)])+20] [times([indexSTOP(1)])+20]], [freqs(freqind)-3 freqs(freqind)-3],'Go', 'fontsize',20);
text([[times([indexSTOP(2)])+20] [times([indexSTOP(2)])+20]], [freqs(freqind)-5 freqs(freqind)-5],'Stop', 'fontsize',30);
text([[times([indexSTOP(3)])+20] [times([indexSTOP(3)])+20]], [freqs(freqind)-7 freqs(freqind)-7],'SSRT', 'fontsize',30);
ylabel('frequency (Hz)','fontsize',30)
xlabel('time (ms)','fontsize',30)
title('successful stop','fontsize',30)
set(gca,'fontsize',30)
 saveas(fu3, [figpath 'STOP_' CLName '_AV_Succ.fig'])
 print(fu3,  '-depsc2 ', [figpath 'STOP_' CLName '_AV_Succ'])

%% plot Fail
 figure ( fu4 );
stat.dimord = 'chan_freq_time';
stat.freq = freqs(freqind1:freqind);
stat.time = times(indexLat(1):indexLat(2));
stat.cfg.channel = {'Chosen'};
stat.label = {'Chosen'}; 
stat.mask = maskedfail;
stat.effect= AVfail;
stat.stat = AVfail;
stat.cfg.equency=[];
j=jet;
colormap(j)
cfg = [];
cfg.channel       = {'Chosen'};
cfg.ylim= [freqs(freqind1) freqs(freqind)];
cfg.zlim           = [-MAXS MAXS];
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
for k = 1:length(STOP_newLat)
line([[times([indexSTOP(k)])] [times([indexSTOP(k)])]],[freqs(freqind1) freqs(freqind)], 'LineWidth',2);
end
line([[times([indexFAIL(3)])] [times([indexFAIL(3)])]],[freqs(freqind1) freqs(freqind)], 'LineStyle', '--', 'LineWidth',2);
text([[times([indexSTOP(1)])+20] [times([indexSTOP(1)])+20]], [freqs(freqind)-3 freqs(freqind)-3],'Go', 'fontsize',20);
text([[times([indexSTOP(2)])+20] [times([indexSTOP(2)])+20]], [freqs(freqind)-5 freqs(freqind)-5],'Stop', 'fontsize',30);
text([[times([indexFAIL(3)])+15] [times([indexFAIL(3)])+15]], [freqs(freqind)-7 freqs(freqind)-7],'RT', 'fontsize',30);
text([[times([indexSTOP(3)])+20] [times([indexSTOP(3)])+20]], [freqs(freqind)-9 freqs(freqind)-9],'SSRT', 'fontsize',30);
ylabel('frequency (Hz)','fontsize',30)
xlabel('time (ms)','fontsize',30)
title('failed stop','fontsize',30)
set(gca,'fontsize',30)
 saveas(fu4, [figpath 'STOP_' CLName '_AV_Failsig.fig'])
 print(fu4,  '-depsc2 ', [figpath 'STOP_' CLName '_AV_Failsig'])

%% plot Go
  figure ( fu5 );
stat.dimord = 'chan_freq_time';
stat.freq = freqs(freqind1:freqind);
stat.time = times(indexLat(1):indexLat(2));
stat.cfg.channel = {'Chosen'};
stat.label = {'Chosen'}; 
stat.mask = maskedgo;
stat.effect= AVgo;
stat.stat = AVgo;
stat.cfg.equency=[];
j=jet;
colormap(j)
cfg = [];
cfg.channel       = {'Chosen'};
cfg.ylim= [freqs(freqind1) freqs(freqind)];
cfg.zlim           = [-MAXS MAXS];
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
for k = 1:length(GO_newLat )
line([[times([indexGo(k)])] [times([indexGo(k)])]],[freqs(freqind1) freqs(freqind)], 'LineWidth',2);
end
text([[times([indexGo(1)])+20] [times([indexGo(1)])+20]], [freqs(freqind)-2 freqs(freqind)-2],'Go', 'fontsize',30);
text([[times([indexGo(2)])+15] [times([indexGo(2)])+15]], [freqs(freqind)-5 freqs(freqind)-5],'RT', 'fontsize',30);
ylabel('frequency (Hz)','fontsize',30)
xlabel('time (ms)','fontsize',30)
title('go','fontsize',30)
set(gca,'fontsize',30)
 saveas(fu5, [figpath 'STOP_' CLName '_AV_Gosig.fig'])
 print(fu5,  '-depsc2 ', [figpath 'STOP_' CLName '_AV_Gosig'])

%% Create difference plots
 
%% subtract common baseline
 
Abs_failT = (ERSPFail(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bfail(:,:,:));
Abs_stopT = (ERSPStop(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bstop(:,:,:));

Abs_novelT = (ERSPNovel(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bnovel(:,:,:));
Abs_standT = (ERSPStand(:,freqind1:freqind,[indexLat(1):indexLat(2)]) + Bstand(:,:,:));

Abs_novel = (Abs_novelT - Bnovelstand_diff);
Abs_stand = (Abs_standT - Bnovelstand_diff);

Abs_fail = (Abs_failT - Bgostop_diff);
Abs_stop = (Abs_stopT - Bgostop_diff);

%% reshape data
Abs_novel = permute(Abs_novel, [2 3 1]);
Abs_stand = permute(Abs_stand, [2 3 1]);

Abs_fail = permute(Abs_fail, [2 3 1]);
Abs_stop = permute(Abs_stop, [2 3 1]);

%% compute significance
NovelStand = {Abs_novel Abs_stand};
FailStop = {Abs_fail Abs_stop};

[stats, df, pNovSt, surrog] = statcond( NovelStand, 'mode', 'bootstrap', 'naccu', 2000);
[stats, df, pFailS, surrog] = statcond( FailStop, 'mode', 'bootstrap', 'naccu', 2000);
 
%% correct for multiple comparisons
[pNoSt,pmaskedNovelStand] = fdr(pNovSt,0.05);
[pFailStop,pmaskedFailStop] = fdr(pFailS,0.05);

Abs_novel = permute(Abs_novel, [3 1 2]);
Abs_stand = permute(Abs_stand, [3 1 2]);

Abs_fail = permute(Abs_fail, [3 1 2]);
Abs_stop = permute(Abs_stop, [3 1 2]);

%% subtract Novel - Standard ERSP matrix
AV_NovelStand = squeeze(mean(Abs_novel-Abs_stand, 1));

%% subtract Stop-Fail ERSP matrix
AV_FailStop = squeeze(mean(Abs_stop - Abs_fail, 1));

%% mask ERSPs
maskNovelStand = (AV_NovelStand.*pmaskedNovelStand);
maskFailStop = (AV_FailStop.*pmaskedFailStop);

%% set maximum for scaling image
a = max(abs([AV_NovelStand]));
MAXNovel  = floor(max(a))+0.2;

b = max(abs([AV_FailStop]));
MAXStop = floor(max(a));

%% reshape data for plotting
AVnovelstand=reshape(AV_NovelStand,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskednovelstand = reshape(pmaskedNovelStand, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

AVfailstop=reshape(AV_FailStop,[1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);
maskedfailstop = reshape(pmaskedFailStop, [1 length(freqs(freqind1:freqind)) length(times(indexLat(1):indexLat(2)))]);

fu6 = figure ( );
fu7 = figure ( );

%% plot Novel - Standard
figure ( fu6 );
stat.dimord = 'chan_freq_time';
stat.freq = freqs(freqind1:freqind);
stat.time = times(indexLat(1):indexLat(2));
stat.cfg.channel = {'Chosen'};
stat.label = {'Chosen'}; 
stat.mask = maskednovelstand ;
stat.effect= AVnovelstand;
stat.stat = AVnovelstand;
stat.cfg.equency=[];
j=jet;
colormap(j)
cfg = [];
cfg.channel       = {'Chosen'};
cfg.ylim= [freqs(freqind1) freqs(freqind)];
cfg.zlim           = [-MAXStop MAXStop];
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
text([[times([indexNOVEL(1)])+15] [times([indexNOVEL(1)])+15]], [freqs(freqind)-2 freqs(freqind)-2],'Sound', 'fontsize',30);
text([[times([indexNOVEL(2)])+15] [times([indexNOVEL(2)])+15]], [freqs(freqind)-3 freqs(freqind)-3],'Target', 'fontsize',30);
ylabel('frequency (Hz)','fontsize',30)
xlabel('time (ms)','fontsize',30)
title('novel - standard','fontsize',25)
set(gca,'fontsize',30)
 saveas(fu6, [figpath 'WM_ ' CLName '_AV_NovelStandDIFF.fig'])
 print(fu6,  '-depsc2 ', [figpath 'WM_' CLName '_AV_NovelStandDIFF'])


%% plot Stop - Fail
 figure ( fu7 );
stat.dimord = 'chan_freq_time';
stat.freq = freqs(freqind1:freqind);
stat.time = times(indexLat(1):indexLat(2));
stat.cfg.channel = {'Chosen'};
stat.label = {'Chosen'}; 
stat.mask = maskedfailstop;
stat.effect= AVfailstop;
stat.stat = AVfailstop;
stat.cfg.equency=[];
j=jet;
colormap(j)
cfg = [];
cfg.channel       = {'Chosen'};
cfg.ylim= [freqs(freqind1) freqs(freqind)];
cfg.zlim           = [-MAXStop MAXStop];
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
for k = 1:length(STOP_newLat)
line([[times([indexSTOP(k)])] [times([indexSTOP(k)])]],[freqs(freqind1) freqs(freqind)], 'LineWidth',2);
end
text([[times([indexSTOP(1)])+20] [times([indexSTOP(1)])+20]], [freqs(freqind)-3 freqs(freqind)-3],'Go', 'fontsize',20);
text([[times([indexSTOP(2)])+15] [times([indexSTOP(2)])+15]], [freqs(freqind)-5 freqs(freqind)-5],'Stop', 'fontsize',30);
text([[times([indexSTOP(3)])+10] [times([indexSTOP(3)])+10]], [freqs(freqind)-7 freqs(freqind)-7],'SSRT', 'fontsize',30);
ylabel('frequency (Hz)','fontsize',30)
xlabel('time (ms)','fontsize',30)
title('successful - failed stop','fontsize',25)
set(gca,'fontsize',30)
 saveas(fu7, [figpath 'STOP_ ' CLName '_AV_StopFailDIFF.fig'])
 print(fu7,  '-depsc2 ', [figpath 'STOP_ ' CLName '_AV_StopFailDIFF'])