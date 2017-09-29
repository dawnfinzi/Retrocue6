%% Script to
% - load the weights and spheres from ICA and assign to the WM and stop
%   files separated
% - compute the dipoles for all ICs (warning: this step takes awhile)
%
% - takes in sbj#_stop_filt_rej_j.set and sbj#_WM_filt_rej_j.set files
% - saves new datasets as S#_STOP_ica.set and S#_WM_ica.set

clear all

[folder, name, ext] = fileparts(which('post_ICA.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

sub = [1 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20]; 

for s = 1:numel(sub)
    
    clear EEG_STOP EEG_WM weights sphere dipoles
    
    % reload individual datasets
    % for each individual dataset reload merged ica weights and spheres and
    % assign to EEG file
    mypath = [folder, '/pp/s', num2str(sub(s)), '/'];
    
    EEG_STOP = pop_loadset(['sbj', num2str(sub(s)), '_stop_filt_rej_j.set'],[mypath]);
    EEG_WM = pop_loadset(['sbj', num2str(sub(s)), '_WM_filt_rej_j.set'],[mypath]);

    load([folder, '/pp/s', num2str(sub(s)), '/ICA', num2str(sub(s)), '.mat']);
    
    % assign weights to stop data
    EEG_STOP.icaweights = weights;
    EEG_STOP.icasphere = sphere;
    EEG_STOP = eeg_checkset(EEG_STOP);
    
    % assign weights to WM data
    EEG_WM.icaweights = weights;
    EEG_WM.icasphere = sphere;
    EEG_WM = eeg_checkset(EEG_WM);
   
    
    %% compute dipoles for all ICs
    EEG_STOP = pop_dipfit_settings( EEG_STOP, 'hdmfile',[eeglabpath, 'plugins/dipfit2.3/standard_BEM/standard_vol.mat'],'coordformat','MNI','mrifile',...
        [eeglabpath, 'plugins/dipfit2.3/standard_BEM/standard_mri.mat'],'chanfile',[eeglabpath, 'plugins/dipfit2.3/standard_BEM/elec/standard_1005.elc',...
        'coord_transform',[0.67469 -15.6307 2.426 0.081417 0.0024349 -1.5728 1.1744 1.0622 1.1487] ,'chansel',[1:EEG_STOP.nbchan] );
    EEG_STOP = pop_multifit(EEG_STOP, [1:size(EEG_STOP.icaweights,1)] ,'threshold',100,'plotopt',{'normlen' 'on'});
    EEG_STOP = eeg_checkset(EEG_STOP);
    
    dipoles = EEG_STOP.dipfit;
    
    EEG_WM.dipfit = dipoles; 
    EEG_WM = eeg_checkset(EEG_WM);
    
    %% save datasets
    EEG_STOP = pop_saveset(EEG_STOP, ['S' num2str(sub(s)) '_STOP_ica.set'], mypath);     
    EEG_WM = pop_saveset(EEG_WM, ['S' num2str(sub(s)) '_WM_ica.set'], mypath); 

end

