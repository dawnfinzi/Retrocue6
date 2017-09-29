%%% script for preprocessing of data, combining of task and then combined
%%% ica (for analysis ala Wessel 2013, 2016 and Wagner under review)

clear all

[folder, name, ext] = fileparts(which('preproc.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

sub = [1 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20]; 

%% process stop data
for s = 1:numel(sub)
    
    savepath = [folder, '/pp/s', num2str(sub(s)), '/'];

    %% load filtered dataset 
    EEG = pop_loadset(['sbj', num2str(sub(s)), '_stop_filt_rej_j.set'],[savepath]);
    
    %% change go markers in SST so that they reflect trial type (SS vs FS)
    if sub(s) == 3
        for i = 1:length(EEG.event)-1
            if EEG.event(i).type == 65532 % if a stop signal trial
                EEG.event(i-1).type = 65500; % recode go signal to distinguish from response first
            end
        end
        % now we can recode based on trial type
        for i = 1:length(EEG.event)-1
            if EEG.event(i).type == 65532 % if a stop signal trial
                if EEG.event(i+1).type == 65531 %if there is a response
                    EEG.event(i-1).type = 65580; % mark go cue as FS
                else
                    EEG.event(i-1).type = 65590; % mark go cue as SS
                end
            end
        end
    else
        for i = 1:length(EEG.event)-1
            if EEG.event(i).type == 65532 % if a stop signal trial
                if EEG.event(i+1).type == 65529 || EEG.event(i+1).type == 65530 %if there is a response
                    EEG.event(i-1).type = 65580; % mark go cue as FS
                else
                    EEG.event(i-1).type = 65590; % mark go cue as SS
                end
            end
        end
    end

    %% go through and reject continuous data by hand
    % save the tmprej file so that you can reload (like below) if neccessary
    load([savepath 'rej_matrix_stop.mat']);

    rejectDataAllABC = TMPREJ;
    [rejectEvents] = eegplot2event(rejectDataAllABC);
    EEG = eeg_eegrej(EEG, rejectEvents);
    EEG = eeg_checkset(EEG);
    eeglab redraw

    % save filtered and artifact rejected dataset
    EEG = pop_saveset(EEG, ['sbj' num2str(sub(s)) '_stop_filt_rej_j.set'], savepath);  
    
end

%% process WM data
for s = 1:numel(sub)
    
    savepath = [folder, '/pp/s', num2str(sub(s)), '/'];

    %% load filtered dataset 
    EEG = pop_loadset(['sbj', num2str(sub(s)), '_WM_filt_rej_j.set'],[savepath]);
    
    %% fix retrocue markers for subject 3
    if sub(s) == 3
        % load behavioral data
        load([folder, '/beh/s3/retrocue_s3.mat']);
        
        trial = 0; 
        for i = 1:length(EEG.event)-1
            if EEG.event(i).type == 65531 % if potentially a tone
                if EEG.event(i-1).type == 65532 % if yes a tone because following a square
                    trial = trial + 1;
                    if data.condition(trial) == 0 % standard
                        EEG.event(i).type = 65530;
                    elseif data.condition(trial) == 1 % novel
                        EEG.event(i).type = 65529;
                    end
                end
            end
        end
    end

    %% save filtered dataset
    EEG = pop_saveset(EEG, ['sbj' num2str(sub(s)) '_WM_filt_j.set'], savepath);  
   
    %% go through and reject continuous data by hand
    % save the tmprej file so that you can reload (like below) if neccessary
    load([savepath 'rej_matrix_WM.mat']);

    rejectDataAllABC = TMPREJ;
    [rejectEvents] = eegplot2event(rejectDataAllABC);
    EEG = eeg_eegrej(EEG, rejectEvents);
    EEG = eeg_checkset(EEG);
    eeglab redraw

    % save filtered and artifact rejected dataset
    EEG = pop_saveset(EEG, ['sbj' num2str(sub(s)) '_WM_filt_rej_j.set'], savepath);  

    
end

%% merge datasets
for s = 1:numel(sub)
    
    savepath = [folder, '/pp/s', num2str(sub(s)), '/'];

     %% merge the two preprocessed datasets
    close all
    eeglab
    EEG = pop_loadset(['sbj', num2str(sub(s)), '_WM_filt_rej_j.set'],[savepath]);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = pop_loadset(['sbj', num2str(sub(s)), '_stop_filt_rej_j.set'],[savepath]);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = pop_mergeset( ALLEEG, [1  2], 0);

    % reject artificially created epochs (from Johanna's script)
    EEG = eeg_regepochs(EEG, 0.5, [0 0.5], NaN); %%epochen von =0.5 sec
    EEG = eeg_checkset(EEG);
    EEG = pop_jointprob(EEG,1,[1:size(EEG.nbchan,1)],4,4,0,1); % berechnet die wahrscheinlichkeit von daten lokaler und globaler wert, Standard deviation 5
    EEG = eeg_checkset(EEG);

    % save as sbj#_merged
    EEG = pop_saveset(EEG, ['sbj' num2str(sub(s)) '_merged.set'], savepath);  
end