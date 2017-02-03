clear all
close all
clc

path = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\';

  subject_folder = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'...
                     'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976'};

  allsubj=   {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'...
                 'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976'};

%%
subject_folder_Neglect = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};
allsubj_Older_Neglect = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};

subject_folder_Control = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976'};
allsubj_Control = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976'};
%%
Monash_Portable_EEG = {'HN018'};%List participants tested using the Bellgrove Lab's portable EEG system

Monash_Stationary_EEG = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976'}; %List participants tested using the Bellgrove Lab's stationary EEG system

%%
No_EyeTracker = {'HN018'};%List participants tested using the Bellgrove Lab's portable EEG system

EyeTracker = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976'}; %List participants tested using the Bellgrove Lab's stationary EEG system

%%

duds = []; %
single_participants = [];



if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
end

% side,motion,ITI
targcodes = zeros(2,2,3); %DN:([side [1=left 2=right] x motion [1=up 2=down] x ITI [inter-target-interval, 3 levels)
targcodes(1,1,:) = [101 105 109]; % left patch, up motion
targcodes(1,2,:) = [102 106 110]; % left patch, down motion
targcodes(2,1,:) = [103 107 111]; % right patch, up motion
targcodes(2,2,:) = [104 108 112]; % right patch, down motion

rtlim=[0.2 2]; %RT must be between 200ms and 2000ms
fs=500; % sample rate %500Hz

BLint = [-100 0];   % baseline interval in ms

ts = -1*fs:2*fs;% in sample points, the ERP epoch
t = ts*1000/fs;

trs = [-.700*fs:fs*.500];% in sample points, the response locked ERP epoch
tr = trs*1000/fs;

master_matrix_R = []; % This saves the matrix for SPSS/R analysis.
total_numtr = 0;

ID_vector=cell(length(subject_folder)*508*length(tr),1); %this will save the subjects ID for each single trial can be pasted into SPSS for ID column. Code at the end of the script clear the emplt cells

HPF=0; %Use high-pass filtered erp? 1=yes, 0=no


current=1;
for s=1:length(allsubj)
    erp_CSDr=[]; erpr=[]; erp_CSD=[]; erp=[]; Pupilr=[]; GAZE_Xr=[]; GAZE_Yr=[]; GAZE_X=[]; GAZE_Y=[];
%     disp(['Subject: ',num2str(s)])
    disp(['Subject: ',allsubj{s}])
    
     if ismember(subject_folder{s},Monash_Stationary_EEG)
        mat_file='_8_to_13Hz_neg1000_to_2000_ARchans1to65_35HzLPF_point0HzHPF_ET.mat';
        ch_N2i = [23;27];
        ch_N2c = [27;23]; % right hemi channels for left target, vice versa.
        ch_for_ipsicon(1,:) = [27;23];
        ch_for_ipsicon(2,:) = [23;27];
        ch_l = [23];
        ch_r = [27];
        ch_front = 5; %5=Fz
        ch_CPP = [25];%25=Pz; 53=CPz
    elseif ismember(subject_folder{s},Monash_Portable_EEG)
        mat_file='_8_to_13Hz_neg1000_to_2000_ARchans1to32_35HzLPF_point0HzHPF_ET.mat';
        ch_N2i = [15;20];
        ch_N2c = [20;15]; % right hemi channels for left target, vice versa.
        ch_for_ipsicon(1,:) = [20;15];
        ch_for_ipsicon(2,:) = [15;20];
        ch_l = [15];
        ch_r = [20];
        ch_front = 2;
        ch_CPP = [13];
    else
        keyboard %participant should be listed in Monash_Portable_EEG or Monash_Stationary_EEG
     end
    
        %Some of the stroke Participants have different focal electrodes for N2c and CPP 
    if strcmp(subject_folder(s),'HN009')
        ch_N2i = [56;27];%56=P5
        ch_N2c = [27;56]; % right hemi channels for left target, vice versa.
        ch_for_ipsicon(1,:) = [27;56];
        ch_for_ipsicon(2,:) = [56;27];
        ch_CPP = [58];%58=P2;
    elseif strcmp(subject_folder(s),'HN019')
        ch_N2i = [56;27];%56=P5
        ch_N2c = [27;56]; % right hemi channels for left target, vice versa.
        ch_for_ipsicon(1,:) = [27;56];
        ch_for_ipsicon(2,:) = [56;27];
        ch_CPP = [20];%20=P5
    end
    
    %     load([path subject_folder{s} '\' allsubj{s} '_ROIs.mat']); %load the participant's individualised Alpha ROI electrodes sensitive to Light
    if ismember(subject_folder{s},subject_folder_Neglect)
        group=1;
        load([path 'Monash_Participants\' subject_folder{s} '\' allsubj{s} mat_file])
    
    elseif ismember(subject_folder{s},subject_folder_Control)
        group=2;
        load([path 'Healthy_Older\' subject_folder{s} '\' allsubj{s} mat_file])
    else
        keyboard
    end
    
    %if the final trial was a miss there will be no RT recorded, just need
    %to add a zero for RT in that case
    if length(allRT)<length(allrespLR)
        allRT(length(allRT):length(allrespLR))=0;
    end
    
    
    if HPF %Use high-pass filtered erp?
        erp=erp_HPF;
    end
    
    allTrials=allTrig; % just renamed this because it makes more sense to me to call it trials
    
    %% calculate the response locked ERPs
    
    erpr = zeros(size(erp,1),length(tr),size(erp,3));
    erp_CSDr = zeros(size(erp_CSD,1),length(tr),size(erp_CSD,3));
    Pupilr = zeros(length(tr),size(Pupil,2));
    GAZE_Xr = zeros(length(tr),size(GAZE_X,2));
    GAZE_Yr = zeros(length(tr),size(GAZE_Y,2));
    validrlock = zeros(1,length(allRT)); % length of RTs.
    for n=1:length(allRT);
        [blah,RTsamp] = min(abs(t_crop*fs/1000-allRT(n))); % get the sample point of the RT.
        if RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t_crop) & allRT(n)>0 % is the RT larger than 1st stim RT point, smaller than last RT point.
            erpr(:,:,n) = erp(:,RTsamp+trs,n);
            erp_CSDr(:,:,n) = erp_CSD(:,RTsamp+trs,n);
            if ismember(subject_folder{s},EyeTracker)
            Pupilr(:,n) = Pupil(RTsamp+trs,n);
            GAZE_Xr(:,n) = GAZE_X(RTsamp+trs,n);
            GAZE_Yr(:,n) = GAZE_Y(RTsamp+trs,n);
            end
            validrlock(n)=1;
        end
    end
  
    
   
    %% DN: master_matrix_R columns:
    %1.Subject number 2.
 
    for trial=1:length(allTrials) % get rid of last trigger?
        total_numtr = total_numtr+1;      
        ID_vector(current:current+(length(tr)-1)) = subject_folder(s);
        %% 1. Subject number:
        master_matrix_R(current:current+(length(tr)-1),1) = s;
        %% 2. Neglect or Control (1 or 2):
        if ismember(subject_folder{s},subject_folder_Neglect)
            master_matrix_R(current:current+(length(tr)-1),2) = 1; % (1=Neglect)
        elseif ismember(subject_folder{s},subject_folder_Control)
            master_matrix_R(current:current+(length(tr)-1),2) = 2; % (2=Control)
        else keyboard
        end
        %% 3. total trial number:
        master_matrix_R(current:current+(length(tr)-1),3) = total_numtr;
        %% 4. inter-subject trial number
        master_matrix_R(current:current+(length(tr)-1),4) = trial;
        %% 5. Pupil Diameter:
        if ismember(subject_folder{s},EyeTracker)
        master_matrix_R(current:current+(length(tr)-1),5)=Pupilr(:,trial);
        else
            master_matrix_R(current:current+(length(tr)-1),5)=0;
        end
        %% 6. CPP:
        master_matrix_R(current:current+(length(tr)-1),6)=mean(erpr(ch_CPP,:, trial),1);
        %% 7. Target Side:
        if ismember(allTrials(trial),targcodes(1,:,:))% any left patch targcode. i.e. left target
            master_matrix_R(current:current+(length(tr)-1),7) = 1;
            TargetSide=1;
        else
            master_matrix_R(current:current+(length(tr)-1),7) = 2; %right target
            TargetSide=2;
        end
        %% 8. N2c:
        master_matrix_R(current:current+(length(tr)-1),8)=erpr(ch_N2c(TargetSide,:),:,trial);
        %% 9. N2i: 
        master_matrix_R(current:current+(length(tr)-1),9)=erpr(ch_N2i(TargetSide,:),:,trial);
        %% 10. Time:
         master_matrix_R(current:current+(length(tr)-1),10)=tr;
        %% 11. CPP CSD transformed
        master_matrix_R(current:current+(length(tr)-1),11)=mean(erp_CSDr(ch_CPP,:, trial),1);
        %% 12. N2c CSD transformed:
        master_matrix_R(current:current+(length(tr)-1),12)=erp_CSDr(ch_N2c(TargetSide,:),:,trial);
        %% 13. N2i CSD transformed: 
        master_matrix_R(current:current+(length(tr)-1),13)=erp_CSDr(ch_N2i(TargetSide,:),:,trial);
        %% 14. GAZE_X:
        if ismember(subject_folder{s},EyeTracker)
         master_matrix_R(current:current+(length(tr)-1),14)=GAZE_Xr(:,trial); 
        else
         master_matrix_R(current:current+(length(tr)-1),14)=0;   
        end
        %% 15. GAZE_Y:
        if ismember(subject_folder{s},EyeTracker)
         master_matrix_R(current:current+(length(tr)-1),15)=GAZE_Yr(:,trial);
        else
         master_matrix_R(current:current+(length(tr)-1),15)=0;
        end
        current=current+length(tr);
    end
end
% find empty cells in ID_vector
emptyCells = cellfun(@isempty,ID_vector);
% remove empty cells
ID_vector(emptyCells) = [];

%Save the data in .csv format to be read into R for inferential stats analysis
 csvwrite (['master_matrix_R_Resp_locked_ERP_Neglect.csv'],master_matrix_R)


