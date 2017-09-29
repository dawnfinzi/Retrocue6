%% Compute Event related spectral perturbations for single ICs in RF

% - load subject and IC index for right frontal cluster
% - load data
% - epoch data relative to tone for retrocue and and relative to go signal for stop signal task
% - get latency of stop signal
% - compute Event related spectral perturbations, timewarp stop signal task trials to stop signal delay
% - save ERSP matrixes for right frontal cluster ICs

clear all
clc

[folder, name, ext] = fileparts(which('computeERSP.m'));
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

Cls_SUBJECT = RFCSJ; % load subject numbers

Cls = RFCcomps; % load IC indices 
   
    
for fa = 1:length(Cls_SUBJECT)
    
    IC = Cls(fa);

    Participant_code = Cls_SUBJECT(fa);
    
    clear EEG_STOP EEG_WM 
    
    %% reload individual datasets with ICA weights and dipoles
    str = char(Participant_code);
    subj_num = str(2:end); 
    if subj_num(1) == '0'
        subj_num = subj_num(2);
    end
    loadpath = [folder, '/pp/s', subj_num, '/'];
    
    EEG_STOP = pop_loadset([strcat(Participant_code, '_STOP_brain_ica.set')],[loadpath]);
    EEG_WM = pop_loadset([strcat(Participant_code, '_WM_brain_ica.set')],[loadpath]);
        
    %% epoch datasets relative to the tone for retrocue
    EEG_S = pop_epoch( EEG_WM, {  '65530'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_S = eeg_checkset(EEG_S);
    
    EEG_N = pop_epoch( EEG_WM, {  '65529'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_N = eeg_checkset(EEG_N);

    %% epoch data relative to go signal for stop signal task
    EEG_SUCCESS = pop_epoch( EEG_STOP, {  '65590'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_SUCCESS = eeg_checkset(EEG_SUCCESS); 
    
    EEG_FAIL = pop_epoch( EEG_STOP, {  '65580'  }, [-1  2.5], 'epochinfo', 'yes');
    EEG_FAIL = eeg_checkset(EEG_FAIL); 
    
    if str2num(subj_num) == 3
        EEG_GO = pop_epoch( EEG_STOP, {  '65535'  }, [-1 2.5], 'epochinfo', 'yes'); 
        EEG_GO = eeg_checkset(EEG_GO);
    else
        EEG_GO = pop_epoch( EEG_STOP, {  '65533'  }, [-1  2.5], 'epochinfo', 'yes');
        EEG_GO = eeg_checkset(EEG_GO); 
    end

    %% get eventlatencies for stop signal task for timewarping to SSD

    STOPons = [10 2000];  
    RESP = [50 2000];

    % go trials
    if str2num(subj_num) == 3
        [RESPgo] = eeg_getepochevent(EEG_GO , {'65531'}, RESP, 'latency'); 
    else
        [RESPgo] = eeg_getepochevent(EEG_GO , {'65529' '65530'}, RESP, 'latency'); 
    end

    % successful stop trials
    [STOPsuc] = eeg_getepochevent(EEG_SUCCESS , '65532', STOPons, 'latency'); 

    % failed stop trials
    [STOPfail] = eeg_getepochevent(EEG_FAIL , '65532', STOPons, 'latency'); 

    % create variables for trials
    offset_GO = zeros(size(1:EEG_GO.trials));
    offset_FAIL = zeros(size(1:EEG_FAIL.trials));
    offset_SUCC = zeros(size(1:EEG_SUCCESS.trials)); 

    GO = [offset_GO' RESPgo'];
    FAIL = [offset_FAIL' STOPfail'];
    SUCC = [offset_SUCC' STOPsuc'];
    
    % set average response time and average SSD based on behavioral data
    GO_newLat = [0 536]; 
    STOP_newLat = [0 252];

    %% Calculate ERSPs
    
    [erspGo, itc, powbaseGo, times, freqs, eboot, pboot, tfdataGo] = newtimef...
        (EEG_GO.icaact(IC,:,:), EEG_GO.pnts, [EEG_GO.xmin EEG_GO.xmax]*1000, EEG_GO.srate, 0,...
        'cycles', [3 0.5], 'tlimits', [-1000 2500], 'freqs', [4 60], 'timesout', 800,...
        'timewarp', [GO], 'timewarpms', [GO_newLat], 'baseline',[-1000 0],...
         'plotitc', 'off', 'plotersp', 'off','padratio', 8);


    [erspFail, itc, powbaseFail, times, freqs, erspboot, itcboot, tfdataFail] = newtimef...
        (EEG_FAIL.icaact(IC,:,:), EEG_FAIL.pnts, [EEG_FAIL.xmin EEG_FAIL.xmax]*1000, EEG_FAIL.srate, 0,...
        'cycles', [3 0.5], 'tlimits', [-1000 2500], 'freqs', [4 60], 'timesout', 800,...
        'timewarp', [FAIL], 'timewarpms', [STOP_newLat], 'baseline',[-1000 0],...
         'plotitc', 'off', 'plotersp', 'off','padratio', 8);


    [erspSucc, itc, powbaseSucc, times, freqs, eboot, pboot, tfdataSucc] = newtimef...
        (EEG_SUCCESS.icaact(IC,:,:), EEG_SUCCESS.pnts, [EEG_SUCCESS.xmin EEG_SUCCESS.xmax]*1000,  EEG_SUCCESS.srate, 0,...
        'cycles', [3 0.5], 'tlimits', [-1000 2500], 'timesout', 800, 'freqs', [4 60],...
        'timewarp', [SUCC], 'timewarpms', [STOP_newLat], 'baseline',[-1000 0],...
        'plotitc', 'off', 'plotersp', 'off','padratio', 8);


    [erspStand, itc, powbaseStand, times, freqs, erspboot, itcboot, tfdataStand] = newtimef...
        (EEG_S.icaact(IC,:,:), EEG_S.pnts, [EEG_S.xmin EEG_S.xmax]*1000, EEG_S.srate, 0,...
        'cycles', [3 0.5], 'tlimits', [-1000 2500], 'freqs', [4 60], 'timesout', 800,...
        'baseline',[-1000 0],...
         'plotitc', 'off', 'plotersp', 'off','padratio', 8);


    [erspNovel, itc, powbaseNovel, times, freqs, eboot, pboot, tfdataNovel] = newtimef...
        (EEG_N.icaact(IC,:,:), EEG_N.pnts, [EEG_N.xmin EEG_N.xmax]*1000,  EEG_N.srate, 0,...
        'cycles', [3 0.5], 'tlimits', [-1000 2500], 'timesout', 800, 'freqs', [4 60],...
        'baseline',[-1000 0],...
        'plotitc', 'off','plotersp', 'off', 'padratio', 8);


    ERSPGo(fa,:,:) = erspGo; 
    ERSPFail(fa,:,:) = erspFail;
    ERSPSucc(fa,:,:) = erspSucc;
    ERSPNovel(fa,:,:) = erspNovel;
    ERSPStand(fa,:,:) = erspStand;

    BASEGo(fa,:) = powbaseGo;
    BASEFail(fa,:) = powbaseFail;
    BASESucc(fa,:) = powbaseSucc;
    BASENovel(fa,:) = powbaseNovel;
    BASEStand(fa,:) = powbaseStand;
    

end

save([mypath 'ERSP_RF_TF.mat'],...
    'ERSPGo', 'BASEGo',...
    'ERSPFail', 'BASEFail',...
    'ERSPSucc', 'BASESucc',...
    'ERSPStand', 'BASEStand',...
    'ERSPNovel', 'BASENovel',...
    'times', 'freqs', '-v7.3');