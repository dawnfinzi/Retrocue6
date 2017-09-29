%% Evaluate clusters and find right frontal cluster

% - load STUDY
% - build preclustring array
% - cluster components
% - plot cluster average scalp maps and dipoles
% - find right frontal cluster (manually identified as cluster 10)
% - plot single IC ERSP images and average cluster ERSP of the right frontal cluster
% - save Subject and IC indices in this cluster

clear all

[folder, name, ext] = fileparts(which('clustering.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

mypath = [folder, '/pp/STUDY/'];

%% load study
[STUDY ALLEEG] = pop_loadstudy('filename', 'RETROCUE_STUDY.study', 'filepath', mypath);

% precluster the data
[STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,...
    {'spec' 'npca' 10 'norm' 1 'weight' 3 'freqrange' [3 30] },...
    {'scalp' 'npca' 10 'norm' 1 'weight' 4 'abso' 1},...
     {'erp' 'npca' 10 'norm' 1 'weight' 3 'timewindow' [0 600] },...
    {'dipoles' 'norm' 1 'weight' 12},...
    {'ersp' 'npca' 10 'norm' 1 'weight' 10 'timewindow' [0 600] 'freqrange' [3 20]},...
    {'finaldim' 'npca' 10});
    
% cluster the data
[STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  16); 

% plot sclpmaps for all clusters
STUDY = std_topoplot(STUDY,ALLEEG,'clusters',[2:17]);
STUDY = std_dipplot(STUDY,ALLEEG,'clusters',[2:17]);

% save clustering results
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename','RETROCUE_STUDY.study','filepath',mypath);

%% Manually load the study in eeglab and identify cluster with right lateralized frontal scalp distribution

rightfront = 10; 

% set parameters for ERSP images
STUDY = pop_erspparams(STUDY, 'topotime',NaN,'topofreq',NaN,'timerange',[-500 1000] ,'freqrange',[4 40], 'subbaseline','on' );
STUDY = pop_statparams(STUDY, 'condstats','on','method','bootstrap');

for ICs = 1:length(STUDY.cluster(rightfront).comps)
STUDY = std_erspplot(STUDY,ALLEEG,'clusters',rightfront, 'comps', ICs );
end

STUDY = std_topoplot(STUDY,ALLEEG,'clusters',rightfront, 'plotsubjects', 'on' );
STUDY = std_topoplot(STUDY,ALLEEG,'clusters',rightfront);
STUDY = std_dipplot(STUDY,ALLEEG,'clusters',rightfront, 'plotsubjects', 'on' );

% plot ERSP image
STUDY = std_erspplot(STUDY,ALLEEG,'clusters',rightfront);

% save subject number and IC index in variable
RFCcomps = STUDY.cluster(rightfront).comps;
RFCSJ = STUDY.subject(STUDY.cluster(rightfront).setinds{1});

save([mypath 'Cluster_RF.mat'], 'RFCcomps', 'RFCSJ')
