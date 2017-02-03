%%% runafew
clear all
close all
clc

%% Monash Neglect participants:
% path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Monash_Participants\';
% subject_folder = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};
% allsubj =        {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};
% 
% allblocks = {[1:15],[1:13],[1:10],[1:10],[1:15],[1:2,4:6],[1:15],[1:20],[1:9]}; % 
% %missing HN009_4.mat
% duds = [];  
% single_participants = [1]; %can put 1 participant, or multiple in here. E.g. [3:8] to process participants 3 to 8
% 
% allbadchans = {[16,32,64]% HN003
%                [16,17,29,32,64]%HN004 
%                [1,2,17,22,36,41,46]% HN009
%                [2,12,46,51]% HN015
%                [2,17,22,42,40,36]%HN020
%                [12,17]%HN019
%                [12,46,40]%HN021
%                [5,9,21,27]}; %HN018

 %% Which EEG system was the participant recorded on? 

% Monash_Portable_EEG = {'HN018'};%List participants tested using the Bellgrove Lab's portable EEG system
% 
% Monash_Stationary_EEG = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021'}; %List participants tested using the Bellgrove Lab's stationary EEG system


%% Older Healthy Participants
% path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Healthy_Older\';
% subject_folder = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976','HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968'};
% allsubj = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976','HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968'};
% 
% allblocks = {[1,3:8],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:3,5:9],[1:9],[3:9],[1,3:9],[1:9],[1:9],[1:9],[2:9],[2:9],[2:9]};
% duds = []; 
% file_start = [1];
% single_participants = [];
% 
%        allbadchans = {[41,17,46,22,1,4,12,33,34,36,37,42] %HN999
%                       [16,17]%HN998
%                       [14,41,17,46,22]%HN996
%                       [12,17,28,41,42,51,52]%HN995
%                       [41,17,46]%HN994
%                       [42]% HN993
%                       [41,17,46,22]% HN992
%                       [12,17,22,41,46,52]% HN990
%                       []% HN989
%                       [22,33,40,41,46] %HN988
%                       [17,41,46]% HN987
%                       [17,42]% HN986
%                       [17,25]% HN985
%                       [41,46,52]% HN983
%                       [4,17,22,41,46,52]% HN982
%                       [9]% HN981
%                       [17,22,51]% HN980
%                       []% HN978
%                       [] % HN977
%                       [17,29,34]% HN976
%                       []% HN975
%                       [17,22,28,41]% HN974
%                       [17,22,36]% HN973
%                       [25,33,36]% HN972
%                       [17,46]% HN971
%                       [12,22]% HN970
%                       []% HN969
%                       []};% HN968
% Monash_Portable_EEG = {};%List participants tested using the Bellgrove Lab's portable EEG system
% Monash_Stationary_EEG = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976','HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968'}; %List participants tested using the Bellgrove Lab's stationary EEG system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Younger Healthy Participants
path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Healthy_Younger\';
subject_folder = {'HN899','HN898','HN897', 'HN896', 'HN895', 'HN894','HN893','HN892', 'HN891', 'HN890','HN889','HN888','HN886','HN885','HN884','HN883','HN882','HN881','HN880','HN879','HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'};
allsubj = {'HN899','HN898','HN897', 'HN896', 'HN895', 'HN894','HN893','HN892', 'HN891', 'HN890','HN889','HN888','HN886','HN885','HN884','HN883','HN882','HN881','HN880','HN879','HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'};

allblocks = {[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[2:9],[1:9],[2:9],[1:9],[1:9],[1:9],[2:9],[1:9],[1:9],[1:9],[1,3:9],[2:9],[2:9]};

duds = []; 
file_start = [];
single_participants = [11];

           allbadchans = {[17,41] %HN899
                      [2,26,50,28]%HN898
                      [1,17,32,52]%HN897
                      [43,51]%HN896
                      []%HN895
                      [17,28,41]% HN894
                      [17,41]% HN893
                      [17,41]% HN892
                      [5,16,17]% HN891
                      [41]% HN890
                      [17,41]%HN889
                      [17]%HN888
                      [17,22]%HN886
                      [42]%HN885
                      [17]% HN884
                      [5,7,12,22,27,36,64]% HN883
                      [12,17,22,41,42]% HN882
                      [4,17,33,37,51,52,55]% HN881
                      [12,42,60,29,23,3,37]% HN880
                      []% HN879
                      [35]%HN878
                      [1,8,28]% HN877
                      [35]% HN876
                      []% HN875
                      [28]% HN874
                      [1,9,20]% HN873
                      [17,22,36,41,46]%HN872  
                      [34]% HN871
                      []% HN870
                      [35]%HN868   
                      [1,33,36]};%HN867
Monash_Portable_EEG = {};%List participants tested using the Bellgrove Lab's portable EEG system
Monash_Stationary_EEG = {'HN899','HN898','HN897', 'HN896', 'HN895', 'HN894','HN893','HN892', 'HN891', 'HN890','HN889','HN888','HN886','HN885','HN884','HN883','HN882','HN881','HN880','HN879','HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'}; %List participants tested using the Bellgrove Lab's stationary EEG system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    allblocks([duds]) = [];
    allbadchans([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder([single_participants]);
    allsubj = allsubj([single_participants]);
    allblocks = allblocks([single_participants]);
    allbadchans = allbadchans([single_participants]);
end

%% CSD

if ~isempty(Monash_Stationary_EEG) 
E = textread('chans65_monash.asc','%s');
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
% MapMontage(M);
[G_monash_stationary,H_monash_stationary] = GetGH(M);
end

if ~isempty(Monash_Portable_EEG)
E = textread('chans32_monash.asc','%s');
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
% MapMontage(M);
[G_monash_portable,H_monash_portable] = GetGH(M);
end

%%

for s=1:length(allsubj)
    disp(allsubj{s})
    
    blocks = allblocks{s};
    badchans = allbadchans{s};
    clear paths files matfiles ET_files ET_matfiles; k=0;
    for n=1:length(blocks)
        k=k+1;
            files{k} = [allsubj{s} '_',num2str(blocks(n)) '.vhdr'];
            paths{k} = [path_temp subject_folder{s} '\'];
            matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' num2str(blocks(n)) '.mat'];
            ET_files{k}=[path_temp 'SamplesAndEvents_combined\' allsubj{s} '_' num2str(blocks(n)) '.asc'];
            ET_matfiles{k} = [path_temp subject_folder{s} '\' allsubj{s} '_' num2str(blocks(n)) '_ET.mat'];
    end            

    
    if ismember(subject_folder{s},Monash_Stationary_EEG)
        G_CSD = G_monash_stationary;
        H_CSD = H_monash_stationary;
        Monash_Stationary_EEG_preprocess %Preprocessing for  Bellgrove Lab's stationary EEG system
    elseif ismember(subject_folder{s},Monash_Portable_EEG)
        G_CSD = G_monash_portable;
        H_CSD = H_monash_portable;
        Monash_Portable_EEG_preprocess %Preprocessing for Bellgrove Lab's portable EEG system
    else
        keyboard %the participant needs to be listed above in either Monash_Stationary_EEG or Monash_Portable_EEG
    end
end
