% getChannelVars.m 
% Computes and saves the channel variances for each block, so these
% can be plotted to look for bad channels using check4badchans.m.
% Variance calc from FFT spectrum so can avoid certain frequencies.
% list ALL BDF files in a study in a cell array of cell arrays such that
% files{s}{n} is the name of the nth file of session (subject) s.
% also list session IDs (subject initials or whatever) in cell array sessionID.
clear all
close all
clc

path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Monash_Participants\'; %Monash PC

subject_folder = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};
allsubj =       {'HN003','HN004','HN009','HN015','HN020','HN019','HN021','HN018'};

                         
Monash_Portable_EEG = {'HN018'};%List participants tested using the Bellgrove Lab's portable EEG system

Monash_Stationary_EEG = {'HN003','HN004','HN009','HN015','HN020','HN019','HN021'}; %List participants tested using the Bellgrove Lab's stationary EEG system

Other_System = {};%List participants tested using the BioSemi system at TCD or QBI

duds = [];
single_participants = [];

file_start = 1;

blocks = {[1:15],[1:13],[1:10],[1:10],[1:15],[1:2,4:6],[1:15],[1:20]};

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

nchan = 32;
for s=file_start:length(subject_folder)
    k=0;
        if ismember(subject_folder{s},Other_System)
            block_temp = blocks{s};
            block_temp(find(ismember(block_temp,0))) = [];
            for n=1:length(block_temp)
                k=k+1;
                files{s}{k} = [path_temp subject_folder{s} '\' allsubj{s},'_',num2str(block_temp(n)) '.bdf'];
            end
        elseif ismember(subject_folder{s},Monash_Portable_EEG) | ismember(subject_folder{s},Monash_Stationary_EEG)
            block_temp = blocks{s};
            block_temp(find(ismember(block_temp,0))) = [];
                for n=1:length(block_temp)
                    k=k+1;
                    paths{s}{k} = [path_temp subject_folder{s} '\'];
                    files{s}{k} = [allsubj{s} '_', num2str(block_temp(n)) '.vhdr'];
                end
            
        end
    matfiles{s} = [path_temp subject_folder{s} '\' allsubj{s} '_chanvars.mat'];
end

% how much of the spectrum to use?
speclims = [0.5 40];  % Limits in Hz

%%%%%%%%%%%%%%%%%%%%% From here it's all standard

h = waitbar(0,'Please wait...');
steps = length(allsubj);
step = file_start-1;                
for s=file_start:length(allsubj)    
    if ismember(subject_folder{s},Other_System)
        chanlocs = readlocs('cap64.loc');
    elseif ismember(subject_folder{s},Monash_Portable_EEG)
        chanlocs = readlocs ('actiCAP32_ThetaPhi.elp','filetype','besa'); %DN for actiCAP 32 from portable system
        paths1 = paths{s};
    elseif ismember(subject_folder{s},Monash_Stationary_EEG)
        chanlocs = readlocs ('actiCAP64_ThetaPhi.elp','filetype','besa'); %DN for actiCAP 64 stationary system
        paths1 = paths{s};
    end
    files1 = files{s};
    matfiles1 = matfiles{s};
    disp(s)
    disp(allsubj{s})
    tic
    if s==file_start
        waitbar(step/steps,h)
    else
        min_time = round((end_time*(steps-step))/60);
        sec_time = round(rem(end_time*(steps-step),60));
        waitbar(step/steps,h,[num2str(min_time),' minutes remaining'])
    end        
    step=step+1;
        
    clear chanVar
    for b=1:length(files1)
        % For the purposes of looking for bad channels, it seems most sensible to leave the BDF referenced as it was recorded.
        % If we average-reference, a bad channel's badness is diluted and may spread to other channels.
        % With a single reference channel, it would be ok, as long as that channel is clean.
        if ismember(subject_folder{s},Other_System)
            EEG = pop_biosig(files1{b}, 'blockepoch', 'off','channels',[1:nchan]);
        elseif ismember(subject_folder{s},Monash_Portable_EEG) | ismember(subject_folder{s},Monash_Stationary_EEG)
            EEG = pop_loadbv(paths1{b},files1{b});
            EEG = pop_rmdat(EEG, {'boundary'},[0 1] ,1); %DN: this deletes the data from 0 to 1 sec around the DC Correction triggers that BrainVision puts in
            for i=1:length([EEG.event.latency]) %DN: round the event latencies back to whole numbers (intergers), because the pop_rmdat line (above) makes non-interger latencies at the points where you kicked out the bad DC Correction data
                EEG.event(i).latency = round([EEG.event(i).latency]);
            end
            EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
        end
        EEG.data = detrend(EEG.data')'; % GL: NB to prevent trending seeming like a bad channel
        
        % Fish out the event triggers and times
        clear trigs stimes
        for i=1:length(EEG.event)
            trigs(i)=EEG.event(i).type;
            stimes(i)=EEG.event(i).latency;
        end
        temp = abs(fft(EEG.data(:,stimes(1):stimes(end))'))'; % FFT amplitude spectrum
        tempF = [0:size(temp,2)-1]*EEG.srate/size(temp,2); % Frequency scale
        chanVar(:,b) = mean(temp(:,find(tempF>speclims(1) & tempF<speclims(2))),2);       % ROW of variances        
    end
    save(matfiles1,'chanVar')
    end_time = toc;
end
close(h)

