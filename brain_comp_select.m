%% Script to help with selecting potential brain components for each subject on eeglab

clear all
close all

[folder, name, ext] = fileparts(which('brain_comp_select.m'));
cd([folder '/'])
home;
current_folder = pwd;
addpath(genpath(current_folder))

mypath = ([current_folder '/']);
eeglabpath = ([current_folder '/eeglab13_6_5b/']);
addpath(eeglabpath)
eeglab

sub = [20];  % change to the subject number that you want to review
s = 1;

mypath = [folder, '/pp/s', num2str(sub(s)), '/'];
    
EEG = pop_loadset(['S' num2str(sub(s)) '_STOP_ica.set'],[mypath]);
 
%% epoch and check components for stop data - to pick brain components manually
EEG = pop_epoch( EEG, {  '65532'  }, [-1  2], 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes'); % epoch around stop signal for selecting comps
eeglab redraw

pop_selectcomps(EEG, [1:size(EEG.icaweights,1)] ); % view components

% view the dipole locations (if RV > 15% - don't include)
pop_dipplot( EEG,[1:size(EEG.icaweights,1)] ,'mri',[eeglabpath,'plugins/dipfit2.3/standard_BEM/standard_mri.mat'],'projlines','on','normlen','on');

% record the numbers of the brain components you want to keep 