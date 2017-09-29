%% Add response error to event codes in retrocue task 
% (this allows for analysis of behavior with eeg later)

clear all

[folder, name, ext] = fileparts(which('add_error_codes.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

sub = [1 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20];  

%%
for s = 1:numel(sub)   
    
    % reload individual datasets
    mypath = [folder, '/pp/s', num2str(sub(s)), '/'];
    behpath = [folder, '/beh/s', num2str(sub(s)), '/'];
    
    EEG = pop_loadset(['sbj', num2str(sub(s)), '_WM_filt_j.set'],[mypath]);
    load([behpath, 'retrocue_s', num2str(sub(s)), '.mat']);
    
    % set up a column for outliers as well
    % (outliers identified in behavioral data as > 2.5 SD from condition
    % mean)
    for i = 1:360
        if ismember(i,data.outliers)
            data.out(i,:) = 1;
        else
            data.out(i,:) = 0;
        end
    end

    %% add error and outlier status to event codes for the tone
    t = 1; 
     if sub(s) == 20
        for i = 1:length(EEG.event)
           if t == 94 % to account for trial 94 missing in eeg data
                t = 95;
           end
           if EEG.event(i).type == 65529 || EEG.event(i).type == 65530 % if a tone
               EEG.event(i).error = data.response_error(t);
               EEG.event(i).outlier = data.out(t);
               t = t + 1; 
           end
        end
     else
        for i = 1:length(EEG.event)
           if EEG.event(i).type == 65529 || EEG.event(i).type == 65530 % if a tone
               EEG.event(i).error = data.response_error(t);
               EEG.event(i).outlier = data.out(t);
               t = t + 1; 
           end
        end
     end

    %% reject WM data
    load([folder, '/pp/s', num2str(sub(s)), '/rej_matrix_WM.mat']);
    rejectDataAllABC = TMPREJ;
    [rejectEvents] = eegplot2event(rejectDataAllABC);
    EEG = eeg_eegrej(EEG, rejectEvents);
    EEG = eeg_checkset(EEG);
    
    %% save datasets
    EEG = pop_saveset(EEG, ['sbj', num2str(sub(s)), '_WM_filt_rej_j.set'], mypath);  
end