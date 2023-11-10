% Script messing with Poisson GLM used by Brooks
% 
% vwadia Feb2023

%% 

setDiskPaths

load([diskPath filesep 'Object_Screening' filesep 'AllITCells_500stim_Scrn.mat'])

m_psths = psths;
m_responses = responses;
m_strctCells = strctCells;

load([diskPath filesep 'Object_Screening' filesep 'AllITCells_500stim_ReScreen.mat'])

a_psths = psths;
a_responses = responses;
a_strctCells = strctCells;

strctCells = [m_strctCells a_strctCells];
psths = [m_psths; a_psths];
responses = [m_responses(:, 1:3); a_responses(:, 1:3)];

ITCellsOnly = 1;

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

timelimits = [-0.17 0.53]; %sec
stimDur = 267;

%% for response latency testing - cells with weird psths that this alg should get 
% 
weirdCells = [2150 4516 949 1056 1383 3762 1087 661 713 1905 1912 782 204 1755 3178 625 651 1036 1981];
allCells = cell2mat(strctCELL(:, 1));

weirdIds = find(allCells == weirdCells);
weirdIds = mod(weirdIds, length(strctCells));

weirdIds(weirdIds == 0) = length(strctCells); 

strctCells = strctCells(weirdIds); 
psths = psths(weirdIds, :);
responses = responses(weirdIds, :);



%% separate cells by session 
% setDiskPaths
% basePath = [diskPath filesep 'Object_Screening' filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115'];
% load([basePath filesep 'PsthandResponses']);
% load([basePath filesep 'strctCells']);
% timelimits = [-0.17 0.53]; %sec
% stimDur = 267;
% 
% psths = screeningPsth;
% 
% strctCELL = struct2cell(strctCells');
% strctCELL = strctCELL';
% 
% % sessID = 'P73CS_ParamObj';
% sessID = 'P82CS_CL_1';
% 
% sessCells = cellfun(@(x) strcmp(x, sessID), strctCELL(:, 8), 'UniformOutput', false);
% 
% sessCells = cell2mat(sessCells);
% 
% strctCells = strctCells(sessCells);
% psths = psths(sessCells, :);
% responses = responses(sessCells, :);



%% Setup orders

% set up orders
stimDir = '500Stimuli';
imageIDs = [1:500]';
order1 = repelem(imageIDs, 6);
order2 = repelem(imageIDs, 4);
catOrd = {};

% afternoon session that got fucked
a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];
a_strct = load([a_basePath filesep 'PsthandResponses']);
order3 = a_strct.order;
clearvars a_strct


% make cat labels 
anovaType = 'CategoryObject';
faceInds = 134:210;
objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
textInds = [264:282 290 400:408];
vegInds = [211:235 357:399];
animInds = [1:84 256:263];

for i = 1:3
    
    if i == 1
        order = order1; % 3000 elements
    elseif i == 2
        order = order2; % 2000 elements
    elseif i == 3
        order = order3;
    end
    
    catOrder = zeros(length(order), 1);
    catOrder(ismember(order, faceInds)) = 1;
    catOrder(ismember(order, textInds)) = 2;
    catOrder(ismember(order, vegInds)) = 3;
    catOrder(ismember(order, animInds)) = 4;
    catOrder(ismember(order, objInds)) = 5;
    
    catOrd{i} = catOrder;
end


%% poisson latency function testing 
% respLat = {};
% for cellIndex = 1:length(strctCells)
%      if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
%         numReps = 6;
%         sortedOrder = order1;
%         labels = catOrd{1};
%         timelimits = [-0.17, 0.33];
%         stimOffDur = 166.6250;
%         stimDur = 166.6250;
%     elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
%         sortedOrder = order3;
%         labels = catOrd{3};
%         timelimits = [-0.17, 0.33];
%         stimOffDur = 133.4680;
%         stimDur = 266.6250;
%     else
%         numReps = 4;
%         sortedOrder = order2;
%         labels = catOrd{2};
%         timelimits = [-0.17, 0.53]; %changed for all cells Oct2022 vwadia
%         stimOffDur = 133.4680;
%         stimDur = 266.6250;
%      end
%     
%      
%     [respLat{cellIndex, 1}, respLat{cellIndex, 2}] = Utilities.computeRespLatPoisson(psths(cellIndex, :), labels, sortedOrder, timelimits, stimDur, true);
%     
%     respLat{cellIndex, 3} = strctCells(cellIndex).Name;
%     
% end
% 
% respCellIds = ~isnan(cell2mat(respLat(:, 1)));
% nonRespCells = strctCells(~respCellIds);
% respCells = strctCells(respCellIds);

%% using p_burst as is

% Open questions:
%     What is anchor? - used to compute cISI which is fed into poisscdf
%     What is the best time window to use for most accurate burst detection?
%         So far (3/10/2023): baseline = stimOn-100 to stimOn same as SFC
%                 period = stimOn+50 to stimOff

% Notes:
%     Significance and min spike # per burst is set inside p_burst
%     Function doesn't handle negative numbers
% tic
 for cellIndex = 1:length(strctCells)
exCell = psths{cellIndex, 1}; % nice animal cell psth

% [BOB, EOB, SOB]=p_burst(InTrain, StartT, StopT,doDisplay,varargin)
% [m1, p1] = max(sum(exCell, 2));
BOB = zeros(1, size(exCell, 1));
EOB = zeros(1, size(exCell, 1));
SOB = zeros(1, size(exCell, 1));
stamps = zeros(size(exCell, 1), 3);
onTimes = [];

% cell baseline firing rate
% avgSpikRate = mean(mean(exCell(:, -timelimits(1)*1e3:-timelimits(1)*1e3+50))); % baseline FR of cell
% avgSpikRate = mean(mean(exCell(:, -timelimits(1)*1e3-100:-timelimits(1)*1e3))); % baseline FR of cell

for it = 1:size(exCell, 1)    
   
    % spike timestamps - note that function gets rid of negative numbers
    times = find(exCell(it, :) == 1);
    
%     if ~isempty(times)
%         keyboard
%     end
    % how does one assess significance? - it's inside the function
    % how does one include a given cell's baseline? - start and stop time are used to compute average FR
    % note start time can't be 0
   
%     [b, e, s] = Utilities.p_burst(times, startT, endT, 0); 

% ----- Current way 09/27/2023
%     can manually input an average spike rate - per trial works best so far
%     avgSpikRate = sum(times > -timelimits(1)*1e3 & times < -timelimits(1)*1e3+50)/50; % baseline FR per trial
%     avgSpikRate = sum(times > -timelimits(1)*1e3-50 & times < -timelimits(1)*1e3)/50; % baseline FR per trial
%      startT = -timelimits(1)*1e3+50;
%     endT = (-timelimits(1)*1e3)+stimDur;
%     [b, e, s] = Utilities.p_burst(times, startT, endT, 0, avgSpikRate); 

% ------ trying different baseline (just of trial period)
    startT = -timelimits(1)*1e3;
    endT = (-timelimits(1)*1e3)+stimDur;
    [b, e, s] = Utilities.p_burst(times, startT, endT, 0); % if you don't feed in average spike rate it will compute it for you

    if ~isempty(b)   
%         train = times(times > -timelimits(1)*1e3);
        train = times;%(times > -timelimits(1)*1e3);
        
        % in cases with multiple burtst take one with max surprise
        if length(s) > 1
            [~, pos] = max(s);
            stamps(it, 1) = train(b(pos));
            stamps(it, 2) = train(e(pos));
            stamps(it, 3) = s(pos);
        else % else take the first 
            stamps(it, 1) = train(b);
            stamps(it, 2) = train(e);
            stamps(it, 3) = s;
        end
    end
end

numRespTrials(cellIndex, 1) = sum(stamps(:, 1) ~= 0);
onTimes = stamps(find(stamps(:, 1) ~= 0), 1);
adjOnTimes = onTimes(onTimes > 170);
adjOnTimes = adjOnTimes - 170;
adj{cellIndex, 1} = adjOnTimes;


end
% toc

respLat(:, 1) = cellfun(@(x) mean(x), adj, 'UniformOutput', false);
respLat(:, 2) = cellfun(@(x) std(x), adj, 'UniformOutput', false);

rL = cell2mat(respLat(:, 1));
rLStd = cell2mat(respLat(:, 2));

%% thresholding after Rlcomp
num_std_dvs = 3.5;
respCellIds = zeros(length(strctCells), 1);
for cellIndex = 1:length(strctCells)
    exCell = psths{cellIndex, 2}; % better for finding threshold crossings
    
    if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
        numReps = 6;
        sortedOrder = order1;
        labels = catOrd{1};
        timelimits = [-0.17, 0.33];
        stimOffDur = 166.6250;
        stimDur = 166.6250;
    elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
        sortedOrder = order3;
        labels = catOrd{3};
        timelimits = [-0.17, 0.33];
        stimOffDur = 133.4680;
        stimDur = 266.6250;
    else
        numReps = 4;
        sortedOrder = order2;
        labels = catOrd{2};
        timelimits = [-0.17, 0.53]; %changed for all cells Oct2022 vwadia
        stimOffDur = 133.4680;
        stimDur = 266.6250;
    end

    if ~isnan(rL(cellIndex))
        m_b_psth(1, 1) = mean(mean(exCell(:, -timelimits(1)*1e3-50:-timelimits(1)*1e3)));
        m_b_psth(1, 2) = std(mean(exCell(:, -timelimits(1)*1e3-50:-timelimits(1)*1e3)));

         ctr = 1;
         m_t_psth = [];
        for lb = unique(labels)'
            m_b_psth(ctr, 1) = mean(mean(exCell(find(labels == lb), -timelimits(1)*1e3-50:-timelimits(1)*1e3)));
            m_b_psth(ctr, 2) = std(mean(exCell(find(labels == lb), -timelimits(1)*1e3-50:-timelimits(1)*1e3)));
            
            if -timelimits(1)*1e3+rL(cellIndex)+ceil(stimDur) > size(exCell, 2)
                t_psth = exCell(find(labels == lb), -timelimits(1)*1e3+floor(rL(cellIndex)):end);
            else
                t_psth = exCell(find(labels == lb), -timelimits(1)*1e3+floor(rL(cellIndex)):-timelimits(1)*1e3+floor(rL(cellIndex))+ceil(stimDur));
            end
           
            
            m_t_psth(ctr, :) = mean(t_psth, 1);
            
            
            ctr = ctr + 1;
        end
        if -timelimits(1)*1e3+rL(cellIndex)+ceil(stimDur) > size(exCell, 2)
            full_t_psth = exCell(:, -timelimits(1)*1e3+floor(rL(cellIndex)):end);
        else
            full_t_psth = exCell(:, -timelimits(1)*1e3+floor(rL(cellIndex)):-timelimits(1)*1e3+floor(rL(cellIndex))+ceil(stimDur));
        end
        
        ctr = 1;
        

        % find max group
        [max_gr, max_val] = Utilities.findMaxGroup(full_t_psth, labels); % find pos group
        [min_gr, min_val] = Utilities.findMaxGroup(full_t_psth, labels, true); % find neg group
                
        if max_val > (m_b_psth(max_gr, 1) + num_std_dvs*m_b_psth(max_gr, 2))...
                ||  min_val < (m_b_psth(min_gr, 1) - num_std_dvs*m_b_psth(min_gr, 2))...
                
                
            respCellIds(cellIndex) = 1; 
        end                
    end
end
 
respCellIds = logical(respCellIds);
nonRespCells = strctCells(~respCellIds);
respCells = strctCells(respCellIds);


%% remove cells with mean FR below 1Hz

for cellIndex = 1:length(strctCells)
    if ~isnan(rL(cellIndex))
        exCell = psths{cellIndex, 2};
        if -timelimits(1)*1e3+ceil(rL(cellIndex))+ceil(stimDur) > size(exCell, 2)
            spike_psth = exCell(:, -timelimits(1)*1e3+floor(rL(cellIndex)):end);
        else
            spike_psth = exCell(:, -timelimits(1)*1e3+floor(rL(cellIndex)):-timelimits(1)*1e3+floor(rL(cellIndex))+ceil(stimDur));
        end
        if mean(mean(spike_psth)) < 0.5
%             disp(cellIndex)
            respCellIds(cellIndex) = 0;
        end
        
    end
    
    
end
    
respCellIds = logical(respCellIds);
nonRespCells = strctCells(~respCellIds);
respCells = strctCells(respCellIds);

%% save those that have v strong responses to a few stim that might get averaged out
for cellIndex = 1:length(strctCells)
    exCell = psths{cellIndex, 1};
    
     if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
        sortedOrder = order1;
        stimDur = 166.6250;
    elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
        sortedOrder = order3;       
        stimDur = 266.6250;
    else
        sortedOrder = order2;
        stimDur = 266.6250;
    end
    
    if ~isnan(rL(cellIndex))
        
        cell_b_psth(1, 1) = mean(mean(exCell(:, -timelimits(1)*1e3-50:-timelimits(1)*1e3)));
        cell_b_psth(1, 2) = std(mean(exCell(:, -timelimits(1)*1e3-50:-timelimits(1)*1e3)));
        ctr = 1;
        notime_m_t_psth = [];
        for sO = unique(sortedOrder)'
            if -timelimits(1)*1e3+rL(cellIndex)+ceil(stimDur) > size(exCell, 2)
                notime_t_psth = exCell(find(sortedOrder == sO), -timelimits(1)*1e3+floor(rL(cellIndex)):end);
            else
                notime_t_psth = exCell(find(sortedOrder == sO), -timelimits(1)*1e3+floor(rL(cellIndex)):-timelimits(1)*1e3+floor(rL(cellIndex))+ceil(stimDur));
            end
            notime_m_t_psth(ctr) = mean(mean(notime_t_psth));
            ctr = ctr+1;
        end
        
        respComp{cellIndex, 1} = (notime_m_t_psth - cell_b_psth(1, 1))./cell_b_psth(1, 2);
    end
end

absmaxtr = cellfun(@(x) max(abs(x)), respComp, 'UniformOutput', false);

% % CRITICAL otherwise the empty cells get removed upon conversion to mat
emptyIds = cellfun(@(x) isempty(x), respComp, 'UniformOutput', false);
absmaxtr(cell2mat(emptyIds)) = {nan};


respCellIds(find(cell2mat(absmaxtr) > 10)) = 1;

respCellIds = logical(respCellIds);
nonRespCells = strctCells(~respCellIds);
respCells = strctCells(respCellIds);


RL = cell2mat(respLat);
RL(~respCellIds, :) = nan;
respLat = num2cell(RL);
%%

% % cells w/ mean FR < 0.5hz in spike collection window
% maybeNRS = strctCells([7 24 33 37 50 152 193 203 206 266 270 291 313 320 321 365 368 373 376 378 380 385 387]);


%% plotting as a visual check
cols = Utilities.distinguishable_colors(length(unique(labels)));
for cellIndex = 1:length(strctCells)
    
    num_std_dvs = 4;
    exCell = psths{cellIndex, 2}; % better for finding threshold crossings
    if exist('weirdCells')
        labels = catOrd{2};
    end
    m_t_psth = [];
    ctr = 1;
    for lb = unique(labels)'
        m_b_psth(ctr, 1) = mean(mean(exCell(find(labels == lb), -timelimits(1)*1e3-100:-timelimits(1)*1e3)));
        m_b_psth(ctr, 2) = std(mean(exCell(find(labels == lb), -timelimits(1)*1e3-100:-timelimits(1)*1e3)));
        t_psth = exCell(find(labels == lb), :);
        %             m_t_psth(ctr, 1) = mean(mean(t_psth));
        m_t_psth(ctr, :) = mean(t_psth, 1);
        ctr = ctr + 1;
    end
    figure;
    hold on;
    title({[num2str(strctCells(cellIndex).ChannelNumber) '\_' num2str(strctCells(cellIndex).Name) '\_' strctCells(cellIndex).brainArea], strctCells(cellIndex).SessionID})
    for i = 1:size(m_t_psth, 1)
        plot(m_t_psth(i, :)', 'Color',  cols(i, :));
        yline(m_b_psth(i, 1) + num_std_dvs*m_b_psth(i, 2), 'Color', cols(i, :));
        yline(m_b_psth(i, 1) - num_std_dvs*m_b_psth(i, 2), 'Color', cols(i, :));
    end
    
    
    xline(-timelimits(1)*1e3-50);
    xline(-timelimits(1)*1e3);
    if ~isnan(rL(cellIndex))
        xline(-timelimits(1)*1e3+rL(cellIndex), 'r')
        xline(-timelimits(1)*1e3+rL(cellIndex)+ceil(stimDur), 'r')
    end
    %         yline(mean(m_t_psth, 2), '--');
    hold off
            pause
    
end

%%
lsidx = find(rLStd == 0);
hsidx = find(rLStd > prctile(rLStd, 99));
nanidx = find(isnan(rLStd));
lowStdCells = strctCells(lsidx);
highStdCells = strctCells(hsidx);
nanCells = strctCells(nanidx);

% test these with threshold - see if one can pull out shit neurons
lowRespTrialCells = strctCells(find(numRespTrials < 4)); 

% restCells = strctCells(setdiff(1:length(strctCells), [lsidx; hsidx; nanidx]'));
% nonRespCells = cat(2, lowStdCells, highStdCells, nanCells);

restCells = strctCells(setdiff(1:length(strctCells), [lsidx; nanidx]'));
nonRespCells = cat(2, lowStdCells, nanCells);
respCells = restCells;


postr = cellfun(@(x) sum(x > 3), respComp, 'UniformOutput', false);
bothtr = cellfun(@(x) sum(abs(x) > 2.5), respComp, 'UniformOutput', false); % 3 is too restrictive here you miss lots of good cells 
maxtr = cellfun(@(x) max(x), respComp, 'UniformOutput', false); % but saved by good score here (8)?
absmaxtr = cellfun(@(x) max(abs(x)), respComp, 'UniformOutput', false);


%% High std, low std and nan std
% outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting' filesep 'nanCells']; maybeNRS = nanCells; 
% outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting' filesep 'lowStdCells']; maybeNRS = lowStdCells; 
% outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting' filesep 'highStdCells']; maybeNRS = highStdCells;
outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting' filesep 'lowFR_trialAvg']; 


if ~exist(outDir, 'dir')
    mkdir(outDir)
end

ims = Utilities.readInFiles([diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'SigandNonSig_Combined'], 'png');
fileNum = nan(length(nonRespCells), 1);
allIds = [1:length(ims)]';

for cIdx = 1:length(maybeNRS)
    
    ftoFind = [maybeNRS(cIdx).brainArea '_' num2str(maybeNRS(cIdx).ChannelNumber) '_' num2str(maybeNRS(cIdx).Name) '_1.png'];
     
    fileNum(cIdx) = structfind(ims, 'name', ftoFind);
    
    ftoMove = [ims(fileNum(cIdx)).folder filesep ims(fileNum(cIdx)).name];
    copyfile(ftoMove, [outDir filesep ims(fileNum(cIdx)).name]);
    
end

%% using p_burst with classwise threshold 
%% Look for classwise threshold crossing across time 

num_std_dvs = 2.5;

% set up orders
stimDir = '500Stimuli';
imageIDs = [1:500]';
order1 = repelem(imageIDs, 6);
order2 = repelem(imageIDs, 4);
catOrd = {};

% afternoon session that got fucked
a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];
a_strct = load([a_basePath filesep 'PsthandResponses']);
order3 = a_strct.order;
clearvars a_strct


% make cat labels 
anovaType = 'CategoryObject';
faceInds = 134:210;
objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
textInds = [264:282 290 400:408];
vegInds = [211:235 357:399];
animInds = [1:84 256:263];

for i = 1:3
    
    if i == 1
        order = order1; % 3000 elements
    elseif i == 2
        order = order2; % 2000 elements
    elseif i == 3
        order = order3;
    end
    
    catOrder = zeros(length(order), 1);
    catOrder(ismember(order, faceInds)) = 1;
    catOrder(ismember(order, textInds)) = 2;
    catOrder(ismember(order, vegInds)) = 3;
    catOrder(ismember(order, animInds)) = 4;
    catOrder(ismember(order, objInds)) = 5;
    
    catOrd{i} = catOrder;
end


respComp = zeros(length(strctCells), 1);

for cellIndex = 1:length(strctCells)
exCell = psths{cellIndex, 1}; % nice animal cell psth

thresholds = [];

if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
    labels = catOrd{1};
    timelimits = [-0.17, 0.33];
    stimOffDur = 166.6250;
    stimDur = 166.6250;
    sortedOrder = order1;
elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
    labels = catOrd{3};
    timelimits = [-0.17, 0.33];
    stimOffDur = 133.4680;
    stimDur = 266.6250;
    sortedOrder = order3;
else
    labels = catOrd{2};
    timelimits = [-0.17, 0.53]; %changed for all cells Oct2022 vwadia
    stimOffDur = 133.4680;
    stimDur = 266.6250;
    sortedOrder = order2;
end

ctr = 1;
threshCross = false;
offset = 50;
binsize = 25;
stepsize = 25;
alpha = 0.01;


b_psth = exCell(:, -timelimits(1)*1e3:-timelimits(1)*1e3+50); % cell baseline Raster
m_b_psth = mean(b_psth, 1);

big_t_psth = exCell(:, -timelimits(1)*1e3+50:-timelimits(1)*1e3+ceil(stimDur));
% windowBegin = -timelimits(1)*1e3+50;
% windowEnd = windowBegin+size(big_t_psth, 2)-binsize;

numBins = length(-timelimits(1)*1e3+stepsize:stepsize:size(big_t_psth, 2)-binsize);
p_sel = nan(numBins, 1);      
lbins = [];
for lb = unique(labels')    
    
%     b_psth = exCell(find(labels == lb), -timelimits(1)*1e3:-timelimits(1)*1e3+50); % cell baseline Raster
%     m_b_psth = mean(b_psth, 2);
    
    t_psth = exCell(find(labels == lb), -timelimits(1)*1e3+50:-timelimits(1)*1e3+ceil(stimDur));
    binCtr = 1;
    for window = 1:stepsize:size(t_psth, 2) % do stepsize = binsize for now to avoid mult comparisons
        
        if window+binsize > size(t_psth, 2)
            schmol_t_psth = t_psth(:, window:end);           
        else
            schmol_t_psth = t_psth(:, window:window+binsize);
        end
        [p_sel(binCtr), ~] = ranksum(m_b_psth, mean(schmol_t_psth, 1));
        
         if p_sel(binCtr) < alpha
             lbins(end+1) = -timelimits(1)*1e3+50+window;
         end
         binCtr = binCtr + 1;

    end
    t = lbins; % list of bins
    
    N = 2; %ceil(binsize/stepsize); % these many consecutive bins have to be significant
    x = diff(t)==stepsize;
    f = find([false,x]~=[x,false]);
    g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
    if ~isempty(t(f(2*g-1)))
        respComp(cellIndex) = 1;
        break
    end

end


% % Should do  test between individual category rasters and cell baseline. 
% multMatrix = Utilities.slidingWindowANOVA(exCell, labels, -timelimits(1)*1e3+offset, alpha, 0, [], binSize, stepsize, 2);
% if ~isempty(multMatrix)
%     
%     % old way
%     %                     respLat = multMatrix{g, 2}; % the first timewindow the groups are significantly different
%     
%     
%     t = cell2mat(multMatrix(:, 2))'; % list of bins
%     
%     N = ceil(25/stepsize); % these many consecutive bins have to be significant
%     x = diff(t)==stepsize;
%     f = find([false,x]~=[x,false]);
%     g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
%     if ~isempty(t(f(2*g-1)))
% %         respLat = multMatrix{f(2*g-1), 2}; % the first timewindow the groups start to be different for a while
%         threshCross = 1; 
%         respComp(cellIndex) = 1;
%     end
%     
%     
% end


% [BOB, EOB, SOB]=p_burst(InTrain, StartT, StopT,doDisplay,varargin)
% [m1, p1] = max(sum(exCell, 2));
% BOB = zeros(1, size(exCell, 1));
% EOB = zeros(1, size(exCell, 1));
% SOB = zeros(1, size(exCell, 1));
% stamps = zeros(size(exCell, 1), 3);
% onTimes = [];


% if threshCross
% 
% for it = 1:size(exCell, 1)
%     
%    
%     % spike timestamps - note that function gets rid of negative numbers
%     times = find(exCell(it, :) == 1);
%     
%     % how does one assess significance? - it's inside the function
%     % how does one include a given cell's baseline? - start and stop time are used to compute average FR
%     % note start time can't be 0
%     startT = -timelimits(1)*1e3+50;
%     endT = (-timelimits(1)*1e3)+stimDur;
% %     [b, e, s] = Utilities.p_burst(times, startT, endT, 0); 
% 
% %     can manually input an average spike rate - per trial works best so far
% %     avgSpikRate = sum(times > -timelimits(1)*1e3 & times < -timelimits(1)*1e3+50)/50; % baseline FR per trial
%     avgSpikRate = sum(times > -timelimits(1)*1e3-100 & times < -timelimits(1)*1e3)/100; % baseline FR per trial
%     [b, e, s] = Utilities.p_burst(times, startT, endT, 0, avgSpikRate); 
%     
%     if ~isempty(b)   
% %         train = times(times > -timelimits(1)*1e3);
%         train = times;%(times > -timelimits(1)*1e3);
%         
%         % in cases with multiple burtst take one with max surprise
%         if length(s) > 1
%             [~, pos] = max(s);
%             stamps(it, 1) = train(b(pos));
%             stamps(it, 2) = train(e(pos));
%             stamps(it, 3) = s(pos);
%         else % else take the first 
%             stamps(it, 1) = train(b);
%             stamps(it, 2) = train(e);
%             stamps(it, 3) = s;
%         end
%     end
% end
% 
% numTrials(cellIndex, 1) = sum(stamps(:, 1) ~= 0);
% onTimes = stamps(find(stamps(:, 1) ~= 0), 1);
% adjOnTimes = onTimes(onTimes > 170);
% adjOnTimes = adjOnTimes - 170;
% adj{cellIndex, 1} = adjOnTimes;
% % train(BOB)
% % train(EOB)
% end

end

% respLat(:, 1) = cellfun(@(x) mean(x), adj, 'UniformOutput', false);
% respLat(:, 2) = cellfun(@(x) std(x), adj, 'UniformOutput', false);

nonRespIds = find(respComp == 0);
respIds = setdiff([1:length(strctCells)], nonRespIds);
nonRespCells = strctCells(nonRespIds); % 48  w/ 20 and 8, 70 w/ 25 and 10, 62 w/ 25 and 8
respCells = strctCells(respIds); % 362 w/ 20 and 8, 340 w/ 25 and 10, 348 w/ 25 and 8


%% using p_burst with individual stim threshold

% things to try
% histogram/scatter of responses to all stimuli - compare max vs min stim and see what the difference is 
% consistency of resplat per stim?? for a visually responsive cell it'll be
% v consistent for the responsive stimuli
% Could also compare to cell baseline not stim baseline

num_std_dvs = 3;

% set up orders
stimDir = '500Stimuli';
imageIDs = [1:500]';
order1 = repelem(imageIDs, 6);
order2 = repelem(imageIDs, 4);
catOrd = {};

% afternoon session that got fucked
a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];
a_strct = load([a_basePath filesep 'PsthandResponses']);
order3 = a_strct.order;
clearvars a_strct



for cellIndex = 1:length(strctCells)
exCell = psths{cellIndex, 1}; % nice animal cell psth

thresholds = [];

if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
    labels = order1;
    timelimits = [-0.17, 0.33];
    stimOffDur = 166.6250;
    stimDur = 166.6250;
elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
    labels = order3;
    timelimits = [-0.17, 0.33];
    stimOffDur = 133.4680;
    stimDur = 266.6250;
else
    labels = order2;
    timelimits = [-0.17, 0.53]; %changed for all cells Oct2022 vwadia
    stimOffDur = 133.4680;
    stimDur = 266.6250;
end

ctr = 1;
threshCross = false;

m_b_psth(1, 1) = mean(mean(exCell(:, -timelimits(1)*1e3:-timelimits(1)*1e3+50)));
m_b_psth(1, 2) = std(mean(exCell(:, -timelimits(1)*1e3:-timelimits(1)*1e3+50)));
% m_b_psth(1, 1) = mean(mean(exCell(:, -timelimits(1)*1e3-50:-timelimits(1)*1e3)));
% m_b_psth(1, 2) = std(mean(exCell(:, -timelimits(1)*1e3-50:-timelimits(1)*1e3)));

% m_b_psth = nan(length(unique(labels)), 2);
m_t_psth = nan(length(unique(labels)), 1);
ctr = 1;
window = 25;
nbins = ceil(stimDur/window);
for lb = unique(labels')   
    
%     b_psth = exCell(find(labels == lb), -timelimits(1)*1e3:-timelimits(1)*1e3+50);
    t_psth = exCell(find(labels == lb), -timelimits(1)*1e3+50:-timelimits(1)*1e3+ceil(stimDur)+50);
    
%     m_b_psth(ctr, 1) = mean(mean(b_psth)); 
%     m_b_psth(ctr, 2) = std(mean(b_psth)); 
%     for nb = 1:nbins
%         if (nb)*window > length(t_psth)
%             schmol_t_psth = t_psth(:, (nb-1)*window+1:end);
%         else
%             schmol_t_psth = t_psth(:, (nb-1)*window+1:nb*window);
%         end
%         m_t_psth(ctr, nb) = mean(mean(t_psth)); 
% 
%     end
    m_t_psth(ctr, 1) = mean(mean(t_psth)); 
    ctr = ctr + 1;

end

respComp{cellIndex} =(m_t_psth - m_b_psth(1, 1))./m_b_psth(1, 2);


end

postr = cellfun(@(x) sum(x > 3), respComp, 'UniformOutput', false)';
bothtr = cellfun(@(x) sum(abs(x) > 2.5), respComp, 'UniformOutput', false)'; % 3 is too restrictive here you miss lots of good cells 
maxtr = cellfun(@(x) max(x), respComp, 'UniformOutput', false)'; % but saved by good score here (8)?
absmaxtr = cellfun(@(x) max(abs(x)), respComp, 'UniformOutput', false)';

% n_respComp = nansum(respComp, 1);
% respLat(:, 1) = cellfun(@(x) mean(x), adj, 'UniformOutput', false);
% respLat(:, 2) = cellfun(@(x) std(x), adj, 'UniformOutput', false);

%%
bt = find(cell2mat(bothtr) < 25);
mt = find(cell2mat(absmaxtr) < 10);

nonRespIds = intersect(bt, mt);
respIds = setdiff([1:length(strctCells)], nonRespIds);
nonRespCells = strctCells(nonRespIds); % 48  w/ 20 and 8, 70 w/ 25 and 10, 62 w/ 25 and 8
respCells = strctCells(respIds); % 362 w/ 20 and 8, 340 w/ 25 and 10, 348 w/ 25 and 8

%% copy to new fodler

% Resp and NonResp
outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting' filesep 'NonRespCells'];
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

ims = Utilities.readInFiles([diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'SigandNonSig_Combined'], 'png');
fileNum = nan(length(nonRespCells), 1);
allIds = [1:length(ims)]';

for cIdx = 1:length(nonRespCells)
    
    ftoFind = [nonRespCells(cIdx).brainArea '_' num2str(nonRespCells(cIdx).ChannelNumber) '_' num2str(nonRespCells(cIdx).Name) '_1.png'];
    
    fileNum(cIdx) = structfind(ims, 'name', ftoFind);
    
    ftoMove = [ims(fileNum(cIdx)).folder filesep ims(fileNum(cIdx)).name];
    copyfile(ftoMove, [outDir filesep ims(fileNum(cIdx)).name]);
    
end

outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting'];

rest = setdiff(allIds, fileNum);
for cIdx = 1:length(rest)
    
    ftoMove = [ims(rest(cIdx)).folder filesep ims(rest(cIdx)).name];

    copyfile(ftoMove, [outDir filesep ims(rest(cIdx)).name]);

end



%% consistency of respLat b/w stim


for cellIndex = 1:length(strctCells)
exCell = psths{cellIndex, 1}; % nice animal cell psth

% [BOB, EOB, SOB]=p_burst(InTrain, StartT, StopT,doDisplay,varargin)
% [m1, p1] = max(sum(exCell, 2));
BOB = zeros(1, size(exCell, 1));
EOB = zeros(1, size(exCell, 1));
SOB = zeros(1, size(exCell, 1));
stamps = zeros(size(exCell, 1), 3);
onTimes = [];

% avgSpikRate = mean(mean(exCell(:, -timelimits(1)*1e3:-timelimits(1)*1e3+50))); % baseline FR of cell
for it = 1:size(exCell, 1)
    
   
    % spike timestamps - note that function gets rid of negative numbers
    times = find(exCell(it, :) == 1);
    
    % how does one assess significance? - it's inside the function
    % how does one include a given cell's baseline? - start and stop time are used to compute average FR
    % note start time can't be 0
    startT = -timelimits(1)*1e3+50;
    endT = (-timelimits(1)*1e3)+stimDur;
%     [b, e, s] = Utilities.p_burst(times, startT, endT, 0); 

%     can manually input an average spike rate - per trial works best so far
%     avgSpikRate = sum(times > -timelimits(1)*1e3 & times < -timelimits(1)*1e3+50)/50; % baseline FR per trial
    avgSpikRate = sum(times > -timelimits(1)*1e3-100 & times < -timelimits(1)*1e3)/100; % baseline FR per trial
    [b, e, s] = Utilities.p_burst(times, startT, endT, 0, avgSpikRate); 
    
    if ~isempty(b)   
%         train = times(times > -timelimits(1)*1e3);
        train = times;%(times > -timelimits(1)*1e3);
        
        % in cases with multiple burtst take one with max surprise
        if length(s) > 1
            [~, pos] = max(s);
            stamps(it, 1) = train(b(pos));
            stamps(it, 2) = train(e(pos));
            stamps(it, 3) = s(pos);
        else % else take the first 
            stamps(it, 1) = train(b);
            stamps(it, 2) = train(e);
            stamps(it, 3) = s;
        end
    end
end


stimStamps = reshape(stamps(:, 1), [4 500]);
meanStim = mean(stimStamps, 1);
stdStim = std(stimStamps);
conAll{cellIndex, 1} = stdStim(meanStim ~= 0);

end

mdCon = cellfun(@(x) mean(x), conAll, 'UniformOutput', false);


