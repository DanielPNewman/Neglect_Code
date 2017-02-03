clear all
close all

%% Monash Neglect participants:
path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Monash_Participants\'; 
subject_folder = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};
allsubj =       {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};
 
blocks = {[1:15],[1:13],[1:10],[1:10],[1:15],[1:2,4:6],[1:15],[1:20]};
 
Monash_Portable_EEG = {'HN018'};%List participants tested using the Bellgrove Lab's portable EEG system

Monash_Stationary_EEG = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021'}; %List participants tested using the Bellgrove Lab's stationary EEG system

Other_System = {};%List participants tested using the BioSemi system at TCD or QBI

%  Younger Healthy Participants
% path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Healthy_Younger\'; 
% subject_folder = {'HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'};
% allsubj = {'HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'}; 

%% Older Healthy Participants
% path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Healthy_Older\';
% subject_folder = {'HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968'};
% allsubj = {'HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968'};

duds = [];

single_participants = [8];

file_start = 1;

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    blocks([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    blocks = blocks(single_participants);
end

for s=file_start:length(subject_folder)
    matfiles{s} = [path_temp subject_folder{s} '\' allsubj{s} '_chanvars.mat'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=file_start:length(subject_folder)
    if ismember(subject_folder{s},Other_System)
        chanlocs = readlocs('cap64.loc');
    elseif ismember(subject_folder{s},Monash_Portable_EEG)
        chanlocs = readlocs ('actiCAP32_ThetaPhi.elp','filetype','besa'); %DN for actiCAP 32 from portable system
    elseif ismember(subject_folder{s},Monash_Stationary_EEG)
        chanlocs = readlocs ('actiCAP64_ThetaPhi.elp','filetype','besa'); %DN for actiCAP 64 stationary system
    end
    disp(subject_folder{s})
    load(matfiles{s})
    % 1,4,7,11
%     chanVar = chanVar(:,1:7);
%     chanVar = chanVar(:,8:14);
%     chanVar = chanVar(:,[1:3,6:14]);
    chanVar = double(chanVar);
    
    badchans = [];  %DN: these two lines will swap and bad channels you identify with the 'changechans' in the next line, just for visualisation
    changechans = []; % must be in same order as badchans.
    chanVar(badchans(1:end),:) = chanVar(changechans(1:end),:);
    
    avVar = mean(chanVar,2); 
        
    % average variance for each channel across all 16 conditions
    % on a second sweep for a given subject, might want to plot topo again
    % after getting rid of a really bad one (to make it easier to see other
    % bad channels) - so do something like:
    % avVar(104) = avVar(103);  % quick hack - make a reall bad chan equal its neighbor
    
    figure;
    topoplot(avVar,chanlocs,'plotchans',[1:length(chanlocs)],'electrodes','numbers');
    title(subject_folder{s})
    
    figure; hold on
    h = plot(chanVar(1:length(chanlocs),:));
    title(subject_folder{s})
    legend(h,'Location','NorthEast');
    
end