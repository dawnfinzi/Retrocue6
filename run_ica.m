%% Run ICA on merged datasets created in combined_notes.m
% recommend not running all subjects in one script (potentially open
% multiple instantiations of matlab and run a few subjects in each) as this
% step is very time consuming 

clear all

[folder, name, ext] = fileparts(which('run_ica.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

sub = [1 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20];  

addpath(eeglabpath)
eeglab

for s = 1:numel(sub)
    mypath = [folder, '/pp/s', num2str(sub(s)), '/'];
    
    bname = (['sbj', num2str(sub(s)), '_merged.set']); 
    EEG = pop_loadset('filename',bname,'filepath',mypath);
    EEG = eeg_checkset( EEG );
    
    [weights, sphere] = binica(EEG.data,'extended',1,'pca',size(EEG.data,1)-1, 'maxsteps', 2000); % run ica
    save([mypath 'ICA' num2str(sub(s)) '.mat'], 'weights', 'sphere')
    clear weights
    clear sphere
end