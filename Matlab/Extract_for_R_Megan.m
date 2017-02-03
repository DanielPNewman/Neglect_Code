clear all
close all
clc

path = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\';
%Healthy Older:
subject_folder = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976','HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968',...
    'HN899','HN898','HN897', 'HN896', 'HN895', 'HN894','HN893','HN892', 'HN891', 'HN890','HN889','HN888','HN886','HN885','HN884','HN883','HN882','HN881','HN880','HN879','HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'};

%Healthy Older:
allsubj=   {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976','HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968',...
    'HN899','HN898','HN897', 'HN896', 'HN895', 'HN894','HN893','HN892', 'HN891', 'HN890','HN889','HN888','HN886','HN885','HN884','HN883','HN882','HN881','HN880','HN879','HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'};


%%
subject_folder_Older = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976','HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968'};
allsubj_Older = {'HN999','HN998','HN996', 'HN995', 'HN994','HN993','HN992','HN990','HN989','HN988','HN987','HN986','HN985','HN983','HN982','HN981','HN980','HN978','HN977','HN976','HN975','HN974','HN973','HN972','HN971','HN970','HN969','HN968'};


subject_folder_Younger = {'HN899','HN898','HN897', 'HN896', 'HN895', 'HN894','HN893','HN892', 'HN891', 'HN890','HN889','HN888','HN886','HN885','HN884','HN883','HN882','HN881','HN880','HN879','HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'};
allsubj_Younger = {'HN899','HN898','HN897', 'HN896', 'HN895', 'HN894','HN893','HN892', 'HN891', 'HN890','HN889','HN888','HN886','HN885','HN884','HN883','HN882','HN881','HN880','HN879','HN878','HN877','HN876','HN875','HN874','HN873','HN872','HN871','HN870','HN868','HN867'};
%%

CSD=0; %Use Current Source Density transformed erp? 1=yes, 0=no

ch_N2i = [23;27];
ch_LR=ch_N2i;
ch_N2c = [27;23]; % right hemi channels for left target, vice versa.
ch_for_ipsicon(1,:) = [27;23];
ch_for_ipsicon(2,:) = [23;27];
ch_l = [23];
ch_r = [27];
ch_front = 5; %5=Fz
ch_CPP = [53];%25=Pz; 53=CPz

fs=500;
ts = -1*fs:2*fs;% in sample points, the ERP epoch
t = ts*1000/fs;

trs = [-.500*fs:fs*.100];% in sample points, the response locked ERP epoch
tr = trs*1000/fs;


% side,motion,ITI
targcodes = zeros(2,2,3); %DN:([side [1=left 2=right] x motion [1=up 2=down] x ITI [inter-target-interval, 3 levels)
targcodes(1,1,:) = [101 105 109]; % left patch, up motion
targcodes(1,2,:) = [102 106 110]; % left patch, down motion
targcodes(2,1,:) = [103 107 111]; % right patch, up motion
targcodes(2,2,:) = [104 108 112]; % right patch, down motion

rtlim=[0.1 2]; %RT must be between 200ms and 2000ms
fs=500; % sample rate %500Hz

master_matrix_R = []; % This saves the matrix for SPSS/R analysis.
total_numtr = 0;
ID_vector=cell(30000,1); %this will save the subjects ID for each single trial can be pasted into SPSS for ID column. Code at the end of the script clear the emplt cells

mat_file='_8_to_13Hz_neg1000_to_2000_ARchans1to65_35HzLPF_point0HzHPF_ET.mat';

for s=1:length(allsubj)
    disp(['Subject: ',num2str(s)])
    disp(['Subject: ',allsubj{s}])
    
    if ismember(subject_folder{s},subject_folder_Older)
        group=1;
        load([path 'Healthy_Older\' subject_folder{s} '\' allsubj{s} mat_file])
        load([path 'Healthy_Older\' subject_folder{s} '\' 'avN2c_ParticipantLevel_peak_amp_index.mat'])
        load([path 'Healthy_Older\' subject_folder{s} '\' 'avN2i_ParticipantLevel_peak_amp_index.mat'])
        load([path 'Healthy_Older\' subject_folder{s} '\' 'ROIs.mat'])
    elseif ismember(subject_folder{s},subject_folder_Younger)
        group=2;
        load([path 'Healthy_Younger\' subject_folder{s} '\' allsubj{s} mat_file])
        load([path 'Healthy_Younger\' subject_folder{s} '\' 'avN2c_ParticipantLevel_peak_amp_index.mat'])
        load([path 'Healthy_Younger\' subject_folder{s} '\' 'avN2i_ParticipantLevel_peak_amp_index.mat'])
        load([path 'Healthy_Younger\' subject_folder{s} '\' 'ROIs.mat'])
    else
        keyboard
    end
    
    LH_ROI=LH_ROI_s;
    RH_ROI=RH_ROI_s;
        
    if CSD
        erp=erp_CSD;
    end
    
    %if the final trial was a miss there will be no RT recorded, just need
    %to add a zero for RT in that case
    if length(allRT)<length(allrespLR)
        allRT(length(allRT):length(allrespLR))=0;
    end
    
    allTrials=allTrig; % just renamed this because it makes more sense to me to call it trials
    
        %% calculate the response locked ERPs
    erpr = zeros(size(erp,1),length(tr),size(erp,3));
    validrlock = zeros(1,length(allRT)); % length of RTs.
    for n=1:length(allRT);
        [blah,RTsamp] = min(abs(t_crop*fs/1000-allRT(n))); % get the sample point of the RT.
        if RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t_crop) & allRT(n)>0 % is the RT larger than 1st stim RT point, smaller than last RT point.
            erpr(:,:,n) = erp(:,RTsamp+trs,n);
            validrlock(n)=1;
        end
    end
    
%%    
    %DN: master_matrix_R columns:
    %Participant(1), ResponseMapping(2), Total Trial no(3),
    %Inter-participant Trial no(4),ITI(5), TargetSide(6)
    %MotionDirection(7), RespLR (8), Accuracy (9), FixBreak/Blink(10), Artefact_pre-target(11),
    %Artefact_post-target(12), Rejected trial (13), Threshold_cohLevel(14) ,
    %RT (15), Pre-target AlphaPower overall (16), Pre-target AlphaPower Left Hemi (17)
    %Pre-target AlphaPower Right Hemi (18), Pre-target AlphaAsym (19),
    %Pre-target Pupil Diameter (20), Number of repeated/invalid trials (21), N2c Amp(22), N2i Amp (23),
	%Post-target Alpha Power Left Hemi (24), Post-target Alpha Power Right Hemi (25),
	%CPP half peak latency (26), N2c peak latency (27), Respose locked CPP slope (28),
	%Older_Younger (29)
    
    for trial=1:length(allTrials) % get rid of last trigger?
        total_numtr = total_numtr+1;
        ID_vector(total_numtr) = allsubj(s);
        %% 1. Subject number:
        master_matrix_R(total_numtr,1) = s;
        %% 2. Response mapping:
            master_matrix_R(total_numtr,2) = 1; %participants always used their right hand to click the left mouse button
        %% 3. total trial number:
        master_matrix_R(total_numtr,3) = total_numtr;
        %% 4. inter-subject trial number
        master_matrix_R(total_numtr,4) = trial;
        %% 5. ITI:
        if ismember (allTrials(trial), targcodes(:,:,1)) % any ITI1 targcode.
            master_matrix_R(total_numtr,5) = 1;
        elseif ismember (allTrials(trial), targcodes(:,:,2))% any ITI2 targcode.
            master_matrix_R(total_numtr,5) = 2;
        else
            master_matrix_R(total_numtr,5) = 3; % any ITI3 targcode.
        end
        %% 6. Target Side:
        if ismember(allTrials(trial),targcodes(1,:,:)) % any left patch targcode. i.e. left target
            TargetSide=1;
            master_matrix_R(total_numtr,6) = 1;
        else
            TargetSide=2;
            master_matrix_R(total_numtr,6) = 2; %right target
        end
        %% 7. Motion Direction:
        if ismember(allTrials(trial),targcodes(:,1,:)) % any up motion targcode.
            master_matrix_R(total_numtr,7) = 1; % Uppward motion
        else
            master_matrix_R(total_numtr,7) = 2; %Downward motion
        end
        %% 8. RespLR - left or right mouse click:
        master_matrix_R(total_numtr,8) = allrespLR(trial); %1=left mouse click; 2=right mouse click; 0= no mouse click
        %% 9. Accuracy:
        if allrespLR(trial)== 0
            master_matrix_R(total_numtr,9) = 0; %no response
        elseif allrespLR(trial)== 1
            master_matrix_R(total_numtr,9) = 1; %correct
        elseif allrespLR(trial)==2  
            master_matrix_R(total_numtr,9) = 2; %they clicked the wrong mouse button
        end
        %% 10. Fixation Break or Blink:
        master_matrix_R(total_numtr,10) = fixation_break_n(trial); %0=no fixation or blink; 1=there was a fixation break or blink
        %% 11. Artefact during pre-target window (-500 - 0ms):
        master_matrix_R(total_numtr,11)= artifact_PretargetToTarget_n(trial);
        %% 12. Artefact during post-target window (-100ms to 100ms after response):
        master_matrix_R(total_numtr,12)= artifact_BLintTo100msPostResponse_n(trial);
        %% 13. Rejected trial
        master_matrix_R(total_numtr,13)=rejected_trial_n(trial);
        %% 14. Threshold_cohLevel 
        %if participant's dot coherence which was determined by the staircase, use this:
%         master_matrix_R(total_numtr,14)=Threshold_cohLevel; 
         %otherwise, set whatever coherence you used (90%)?
        master_matrix_R(total_numtr,14)=90;
        %% 15. Reaction time (RT):
        master_matrix_R(total_numtr,15)=allRT(trial)*1000/fs;
        %% 16. Pre-target Alpha Power overall (combining the two ROIs):
        master_matrix_R(total_numtr,16)=squeeze(mean(mean(Alpha([LH_ROI_s RH_ROI_s],find(Alpha_smooth_time==-500):find(Alpha_smooth_time==0),trial),1),2));
        %% 17. Pre-target Alpha Power Left Hemi:
        master_matrix_R(total_numtr,17)=squeeze(mean(mean(Alpha([LH_ROI_s],find(Alpha_smooth_time==-500):find(Alpha_smooth_time==0),trial),1),2));
        %% 18. Pre-target Alpha Power Right Hemi:
        master_matrix_R(total_numtr,18)=squeeze(mean(mean(Alpha([RH_ROI_s],find(Alpha_smooth_time==-500):find(Alpha_smooth_time==0),trial),1),2));
        %% 19.  Pre-target AlphaAsym:
        master_matrix_R(total_numtr,19)=(master_matrix_R(total_numtr,18)-master_matrix_R(total_numtr,17))/(master_matrix_R(total_numtr,18)+master_matrix_R(total_numtr,17)); %(RightHemiROI - LeftHemiROI)/(RightHemiROI + LeftHemiROI)
        %% 20. Pre-target Pupil Diameter:
        master_matrix_R(total_numtr,20)=mean(Pupil(find(t_crop==-500):find(t_crop==0),trial));
        %% 21. Number of repeated trials due to fixation breaks:
        master_matrix_R(total_numtr,21)=length(allRT(~(~fixation_break_n) & ~(~rejected_trial_n)));  
        %% 22. N2c Amp (using PARTICIPANT LEVEL AVERAGE to define N2c measurement window):
        window=25; %this is the time (in samples) each side of the peak latency - so 25 is actually a 100ms window (since fs=500 and this is done each side of the peak latency)
        master_matrix_R(total_numtr,22)=mean(mean(erp(ch_N2c(TargetSide,:),(avN2c_ParticipantLevel_peak_amp_index_s(TargetSide))-window:avN2c_ParticipantLevel_peak_amp_index_s(TargetSide)+window,trial),1));
        %% 23. N2i Amp (using PARTICIPANT LEVEL AVERAGE to define N2i measurement window):
        master_matrix_R(total_numtr,23)=mean(mean(erp(ch_LR(TargetSide,:),avN2c_ParticipantLevel_peak_amp_index_s(TargetSide)-window:avN2c_ParticipantLevel_peak_amp_index_s(TargetSide)+window,trial),1));
        %% 24. Post-target Alpha Power Left Hemi:
        master_matrix_R(total_numtr,24)=squeeze(mean(mean(Alpha([LH_ROI],find(Alpha_smooth_time==300):find(Alpha_smooth_time==1000),trial),1),2));
        %% 25. Post-target Alpha Power Right Hemi:
        master_matrix_R(total_numtr,25)=squeeze(mean(mean(Alpha([RH_ROI],find(Alpha_smooth_time==300):find(Alpha_smooth_time==1000),trial),1),2));
        %% 26. CPP half peak latency:
        half_max_peak=max(erp(ch_CPP,find(t_crop==0):find(t_crop==1500),trial))/2;

        half_max_peak_index=find(erp(ch_CPP,find(t_crop==0):find(t_crop==1500),trial)>=half_max_peak,1,'first')+length(find(erp(t_crop<0)));
        if half_max_peak<0
            master_matrix_R(total_numtr,26)=0;
        elseif isempty(half_max_peak_index)
                master_matrix_R(total_numtr,26)=0;
        else 
        master_matrix_R(total_numtr,26)=t_crop(half_max_peak_index);
        end
        
        %% 27. N2c peak latency:
        % NB: this is quarter the window size in samples, each sample = 2ms and
        % it's this either side.
        window_size = 25;
        % contra and search frames encompass the entire negativity.
        contra_peak_t = [150,500]; contra_peak_ts(1) = find(t_crop==contra_peak_t(1)); contra_peak_ts(2) = find(t_crop==contra_peak_t(2));
        % SEARCHING FOR PEAK LATENCY
        clear N2c_peak_latencies
        % for each trial...
        % search your search timeframe, defined above by contra_peak ts, in sliding windows
        % NB this is done in samples, not time, it's later converted.
        clear win_mean win_mean_inds
        counter2 = 1;
        for j = contra_peak_ts(1)+window_size:contra_peak_ts(2)-window_size
            % get average amplitude of sliding window from N2pc electrode
            if master_matrix_R(total_numtr,6) == 1; %if left target, measure from right hemi (ch_LR(2,:)) electrodes
                win_mean(counter2) = squeeze(mean(mean(erp(ch_LR(2,:),j-window_size:j+window_size,trial),1),2));
            elseif master_matrix_R(total_numtr,6) == 2; %if right target, measure from left hemi (ch_LR(1,:)) electrodes
                win_mean(counter2) = squeeze(mean(mean(erp(ch_LR(1,:),j-window_size:j+window_size,trial),1),2));
            end
            % get the middle sample point of that window
            win_mean_inds(counter2) = j;
            counter2 = counter2+1;
        end
        % find the most negative amplitude in the resulting windows
        [~,ind_temp] = min(win_mean);
        % get the sample point which had that negative amplitude
        N2pc_min_ind = win_mean_inds(ind_temp);
        
        % if the peak latency is at the very start or end of the search
        % timeframe, it will probably be bogus. set to NaN.
        if ind_temp==1 | ind_temp==length(win_mean)
            master_matrix_R(total_numtr,27)= 0; %%DN: make it 0 instead of NaN, will remove these in R
        else
            % it's good! add it in.
            master_matrix_R(total_numtr,27)=t_crop(N2pc_min_ind);  %N2c_peak_latencies(trial)= t(N2pc_min_ind);
        end
       %% 28. Respose locked CPP slope: (just fitting a straight line, like in Kelly and O'Connel J.Neuro, but on trial-by-trial basis)
        slope_timeframe = [-150,-10];
        if validrlock(trial)
            coef = polyfit(tr(tr>slope_timeframe(1) & tr<slope_timeframe(2)),(erpr(ch_CPP,tr>slope_timeframe(1) & tr<slope_timeframe(2),trial)),1); % coef returns 2 coefficients fitting r = slope * x + intercept
            master_matrix_R(total_numtr,28) = coef(1); %slope
            %                 r = coef(1) .* tr(tr>slope_timeframe(1) & tr<slope_timeframe(2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            %                 figure
            %                 plot(tr,erpr(ch_CPP,:,trial),'color','k');
            %                 hold on;
            %                 plot(tr(tr>slope_timeframe(1) & tr<slope_timeframe(2)), r, ':');
            %                 line(xlim,[0,0],'Color','k');
            %                 line([0,0],ylim,'Color','k');
            %                 line([slope_timeframe(1),slope_timeframe(1)],ylim,'linestyle',':');
            %                 line([slope_timeframe(2),slope_timeframe(2)],ylim,'linestyle',':');
            %                 hold off;
        end
      %% 29. Group (Older (1) vs. Younger(2)) 
       master_matrix_R(total_numtr,29) = group;
      
    end
end
  
% find empty cells in ID_vector
emptyCells = cellfun(@isempty,ID_vector);
% remove empty cells
ID_vector(emptyCells) = [];

%Save the data in .csv format to be read into R for inferential stats analysis
csvwrite ('master_matrix_R_Older_Younger.csv',master_matrix_R)
cell2csv ('ID_vector_Older_Younger.csv',ID_vector)
