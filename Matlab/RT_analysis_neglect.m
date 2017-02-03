%% Coded by Ger Loughnane 12/02/2016
clear
close all
clc

%% Load subjects
path_temp = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\BL_Behav_Data_Files\';

subject_folder = {'HNBL006'};

session_folder = {{'HNBL006_2','HNBL006_3H','HNBL006_4L','HNBL006_5L','HNBL006_6H'},...
    };
session_tags = {{'','H','L','L','H'}};
blocks = {{[1:10],[1:10],[1:10],[11:20],[11:20]}};

allsubj = {'HNBL006'};

duds = [];
single_participants = [];

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    session_folder([duds]) = [];
    session_tags([duds]) = [];
    blocks([duds]) = [];
    allsubj([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    session_folder = session_folder(single_participants);
    session_tags = session_tags(single_participants);
    blocks = blocks(single_participants);
    allsubj = allsubj(single_participants);
end
filenames={}; session_number=[]; block_number=[]; light_condition=[];
for s = 1:length(allsubj)
    k=0;
    for session = 1:length(session_folder{s})
        for n = 1:length(blocks{s}{session})
            k=k+1;
            filenames{s,k} = [path_temp subject_folder{s} '\' ...
                session_folder{s}{session} '\' allsubj{s} '_' session_tags{s}{session} num2str(blocks{s}{session}(n)) '.mat'];
            session_number(s,k) = session;
            block_number(s,k) = n;
            if strcmp(session_tags{s}{session},'')
                light_condition(s,k) = 1;
            elseif strcmp(session_tags{s}{session},'L')
                light_condition(s,k) = 2;
            elseif strcmp(session_tags{s}{session},'H')
                light_condition(s,k) = 3;
            end
        end
    end
end

%% Conditions
% for i = 1:length(conditiondescrip)
%     disp(conditiondescrip{i})
% end
% Trigger 1: coherence 90, patch 1, motion dir 90, coh motion onset 1800
% Trigger 2: coherence 90, patch 1, motion dir 270, coh motion onset 1800
% Trigger 3: coherence 90, patch 2, motion dir 90, coh motion onset 1800
% Trigger 4: coherence 90, patch 2, motion dir 270, coh motion onset 1800
% Trigger 5: coherence 90, patch 1, motion dir 90, coh motion onset 2800
% Trigger 6: coherence 90, patch 1, motion dir 270, coh motion onset 2800
% Trigger 7: coherence 90, patch 2, motion dir 90, coh motion onset 2800
% Trigger 8: coherence 90, patch 2, motion dir 270, coh motion onset 2800
% Trigger 9: coherence 90, patch 1, motion dir 90, coh motion onset 3800
% Trigger 10: coherence 90, patch 1, motion dir 270, coh motion onset 3800
% Trigger 11: coherence 90, patch 2, motion dir 90, coh motion onset 3800
% Trigger 12: coherence 90, patch 2, motion dir 270, coh motion onset 3800
% side,motion,ITI
targcodes = zeros(2,2,3);
targcodes(1,1,:) = [101 105 109]; % left patch, up motion, ITI1:ITI3
targcodes(1,2,:) = [102 106 110]; % left patch, down motion, ITI1:ITI3
targcodes(2,1,:) = [103 107 111]; % right patch, up motion, ITI1:ITI3
targcodes(2,2,:) = [104 108 112]; % right patch, down motion, ITI1:ITI3

rtlims = [0.2 1]; % RT limits for plots below, but ALL RTs go into main matrix.
%% Go through files
for s = 1:length(allsubj)
    numtr=0; allRTs=[]; allPerf=[]; allTrig=[]; allHemi=[]; allMotion=[]; allITI=[]; allSession=[]; allBlockNum=[]; allLightConds =[]; 
    for f = 1:size(filenames,2)
        load(filenames{s,f});
        numcond = length(conditiondescrip);
        cohmo_trigs = find(PTBtrig>100 & PTBtrig<=100+numcond); 
        for n=1:length(cohmo_trigs);
            numtr=numtr+1;
            % RTs
            allRTs(numtr) = 0;
            stimtime = PTBtrigT(cohmo_trigs(n));
            %   nextresp = find(RespT>stimtime,1);
            nextresp = find(RespT>stimtime & RespT<stimtime+3,1); %DN: just put a limit (3000ms) on  response times to be included so any response time over 2500ms is kicked out for this analysis
            if ~isempty(nextresp)
                allRTs(numtr) = RespT(nextresp) - stimtime;
                allPerf(numtr) = 1;
            else
                allPerf(numtr) = 0;
            end
            % Triggers etc.
            allTrig(numtr) = PTBtrig(cohmo_trigs(n));
            if ismember(PTBtrig(cohmo_trigs(n)),targcodes(1,:,:)), allHemi(numtr) = 1;
            elseif ismember(PTBtrig(cohmo_trigs(n)),targcodes(2,:,:)), allHemi(numtr) = 2; 
            end
            if ismember(PTBtrig(cohmo_trigs(n)),targcodes(:,1,:)), allMotion(numtr) = 1;
            elseif ismember(PTBtrig(cohmo_trigs(n)),targcodes(:,2,:)), allMotion(numtr) = 2; 
            end
            if ismember(PTBtrig(cohmo_trigs(n)),targcodes(:,:,1)), allITI(numtr) = 1;
            elseif ismember(PTBtrig(cohmo_trigs(n)),targcodes(:,:,2)), allITI(numtr) = 2;
            elseif ismember(PTBtrig(cohmo_trigs(n)),targcodes(:,:,3)), allITI(numtr) = 3; 
            end
            
            allSession(numtr) = session_number(s,f);
            allBlockNum(numtr) = block_number(s,f);
            allLightConds(numtr) = light_condition(s,f);
        end
    end
    % Data which can be read into SPSS in the following order
    all_data = [allRTs',allPerf',allHemi',allMotion',allITI',allLightConds',allSession',allBlockNum'];
    all_subject_data = [s*ones(length(allRTs),1),allRTs',allPerf',allHemi',allMotion',allITI',allLightConds',allSession',allBlockNum'];
    
%     a = sort(allRTs); figure, plot(a)
    [N1,edges1] = histcounts(allRTs(find(allPerf==1 & allHemi==1 & allRTs>rtlims(1) & allRTs<rtlims(2))),10); 
    for i = 1:length(edges1)-1
        values1(i) = mean(edges1(i:i+1));
    end
    [N2,edges2] = histcounts(allRTs(find(allPerf==1 & allHemi==2 & allRTs>rtlims(1) & allRTs<rtlims(2))),10);
    for i = 1:length(edges2)-1
        values2(i) = mean(edges2(i:i+1));
    end
    clear h
    figure
    h(1) = plot(values1,N1,'LineWidth',2); hold on
    h(2) = plot(values2,N2,'LineWidth',2);
    legend(h,{'Left','Right'})
    title([allsubj{s},': Left RT median = ',num2str(median(allRTs(find(allPerf==1 & allHemi==1)))),...
        ', Right RT median = ',num2str(median(allRTs(find(allPerf==1 & allHemi==2))))])
    xlabel('RT (seconds)'), ylabel('Trial Count')
    figure
    histogram(allRTs)
    xlabel('RT (seconds)'), ylabel('Trial Count')
    title([allsubj{s}])
    
    save([path_temp,subject_folder{s},'\',allsubj{s},'_RTs.mat'],'allRTs','allPerf','allTrig','allHemi','allMotion','allITI','allSession','allBlockNum','allLightConds','all_data');
end
save([path_temp,'\ALL_SUBJ_RTs.mat'],'all_subject_data');

return
%% quick analysis by light condition: NB only works on one subject at a time for now!!! can easily change for several subjects.
for lightcond = 1:3
    for side = 1:2
        RT_mean(side,lightcond) = mean(allRTs(find(allPerf==1 & allHemi==side & allRTs>rtlims(1) & allRTs<rtlims(2) & allLightConds==lightcond)));
    end
    % left minus right/mean(left and right): more negative means faster to
    % the left.
    RT_index(lightcond) = (RT_mean(1,lightcond)-RT_mean(2,lightcond))/((RT_mean(1,lightcond)+RT_mean(2,lightcond))/2);
end
figure
for side = 1:2
    h(side) = plot(RT_mean(side,:),'LineWidth',2); hold on
end
set(gca,'xtick',[1,2,3],'xticklabel',{'No Light','Low Light','High Light'})
xlabel('Light Condition'), ylabel('RT(sec)')
legend(h,{'Left','Right'})
title('RTs by target side and light condition')

figure
plot(RT_index,'LineWidth',2); hold on

set(gca,'xtick',[1,2,3],'xticklabel',{'No Light','Low Light','High Light'},'ylim',[0 max(RT_index)])
xlabel('Light Condition'), ylabel('RT index')
title('RT index by light condition (more positive = more rightward biased)')

%% analysis by session order
for sesh = 1:5
    for side = 1:2
        RT_mean2(side,sesh) = mean(allRTs(find(allPerf==1 & allHemi==side & allRTs>rtlims(1) & allRTs<rtlims(2) & allSession==sesh)));
    end
    % left minus right/mean(left and right): more negative means faster to
    % the left.
    RT_index2(sesh) = (RT_mean2(1,sesh)-RT_mean2(2,sesh))/((RT_mean2(1,sesh)+RT_mean2(2,sesh))/2);
end
figure
for side = 1:2
    h(side) = plot(RT_mean2(side,:),'LineWidth',2); hold on
end
set(gca,'xtick',[1,2,3,4,5],'xticklabel',{'Session 1','Session 2','Session 3','Session 4','Session 5'})
xlabel('Session Number'), ylabel('RT(sec)')
legend(h,{'Left','Right'})
title('RTs by target side and session order')

figure
plot(RT_index2,'LineWidth',2); hold on
set(gca,'xtick',[1,2,3,4,5],'xticklabel',{'Session 1','Session 2','Session 3','Session 4','Session 5'},'ylim',[0 max(RT_index2)])
xlabel('Session Number'), ylabel('RT index')
title('RT index by session order (more positive = more rightward biased)')
