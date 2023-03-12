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


%% for response latency testing - cells with weird psths that this alg should get 

weirdCells = [2150 4516 949 1056 1383 661 713 1905 204 1755 3178 625 651 1036 1981 ];
allCells = cell2mat(strctCELL(:, 1));

weirdIds = find(allCells == weirdCells);
weirdIds = mod(weirdIds, length(strctCells));

weirdIds(weirdIds == 0) = length(strctCells); 

strctCells = strctCells(weirdIds); 
psths = psths(weirdIds, :);
responses = responses(weirdIds, :);

%% using p_burst as is

% Open questions:
%     What is anchor? - used to compute cISI which is fed into poisscdf
%     What is the best time window to use for most accurate burst detection?
%         So far (3/10/2023): baseline = stimOn-100 to stimOn same as SFC
%                 period = stimOn+50 to stimOff

% Notes:
%     Significance and min spike # per burst is set inside p_burst
%     Function doesn't handle negative numbers

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

numTrials(cellIndex, 1) = sum(stamps(:, 1) ~= 0);
onTimes = stamps(find(stamps(:, 1) ~= 0), 1);
adjOnTimes = onTimes(onTimes > 170);
adjOnTimes = adjOnTimes - 170;
adj{cellIndex, 1} = adjOnTimes;
% train(BOB)
% train(EOB)
end

respLat(:, 1) = cellfun(@(x) mean(x), adj, 'UniformOutput', false);
respLat(:, 2) = cellfun(@(x) std(x), adj, 'UniformOutput', false);


%% using p_burst with classwise threshold 

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
for lb = unique(labels')    
    b_psth = exCell(find(labels == lb), -timelimits(1)*1e3:-timelimits(1)*1e3+50);
    t_psth = exCell(find(labels == lb), -timelimits(1)*1e3+50:-timelimits(1)*1e3+ceil(stimDur));
    
    if mean(mean(t_psth)) >= mean(mean(b_psth)) + num_std_dvs*std(mean(b_psth))
        threshCross = true;
    end
    ctr = ctr+1;
end

% [BOB, EOB, SOB]=p_burst(InTrain, StartT, StopT,doDisplay,varargin)
% [m1, p1] = max(sum(exCell, 2));
BOB = zeros(1, size(exCell, 1));
EOB = zeros(1, size(exCell, 1));
SOB = zeros(1, size(exCell, 1));
stamps = zeros(size(exCell, 1), 3);
onTimes = [];


if threshCross

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

numTrials(cellIndex, 1) = sum(stamps(:, 1) ~= 0);
onTimes = stamps(find(stamps(:, 1) ~= 0), 1);
adjOnTimes = onTimes(onTimes > 170);
adjOnTimes = adjOnTimes - 170;
adj{cellIndex, 1} = adjOnTimes;
% train(BOB)
% train(EOB)
end

end

respLat(:, 1) = cellfun(@(x) mean(x), adj, 'UniformOutput', false);
respLat(:, 2) = cellfun(@(x) std(x), adj, 'UniformOutput', false);



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

m_b_psth(1, 1) = mean(mean(exCell(:, -timelimits(1)*1e3-100:-timelimits(1)*1e3)));
m_b_psth(1, 2) = std(mean(exCell(:, -timelimits(1)*1e3-100:-timelimits(1)*1e3)));

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
bothtr = cellfun(@(x) sum(abs(x) > 2.5), respComp, 'UniformOutput', false)'; % score cutoff 20?
maxtr = cellfun(@(x) max(x), respComp, 'UniformOutput', false)'; % but saved by good score here (8)?
absmaxtr = cellfun(@(x) max(abs(x)), respComp, 'UniformOutput', false)';

% n_respComp = nansum(respComp, 1);
% respLat(:, 1) = cellfun(@(x) mean(x), adj, 'UniformOutput', false);
% respLat(:, 2) = cellfun(@(x) std(x), adj, 'UniformOutput', false);
bt = find(cell2mat(bothtr) < 20);
mt = find(cell2mat(absmaxtr) < 8);

nonRespIds = intersect(bt, mt);
respIds = setdiff([1:length(strctCells)], nonRespIds);
nonRespCells = strctCells(nonRespIds);
respCells = strctCells(respIds);

%% copy to new fodler

% outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting' filesep 'NonRespCells'];
% if ~exist(outDir, 'dir')
%     mkdir(outDir)
% end
% 
% ims = Utilities.readInFiles([diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'SigandNonSig_Combined'], 'png');
% fileNum = nan(length(nonRespCells), 1);
% allIds = [1:length(ims)]';
% 
% for cIdx = 1:length(nonRespCells)
%     
%     ftoFind = [nonRespCells(cIdx).brainArea '_' num2str(nonRespCells(cIdx).ChannelNumber) '_' num2str(nonRespCells(cIdx).Name) '_1.png'];
%     
%     fileNum(cIdx) = structfind(ims, 'name', ftoFind);
%     
%     ftoMove = [ims(fileNum(cIdx)).folder filesep ims(fileNum(cIdx)).name];
%     copyfile(ftoMove, [outDir filesep ims(fileNum(cIdx)).name]);
%     
% end
% 
% outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting'];
% 
% rest = setdiff(allIds, fileNum);
% for cIdx = 1:length(rest)
%     
%     ftoMove = [ims(rest(cIdx)).folder filesep ims(rest(cIdx)).name];
% 
%     copyfile(ftoMove, [outDir filesep ims(rest(cIdx)).name]);
% 
% end
outDir = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'AllCells' filesep 'RLTesting' filesep 'Maybe'];
if ~exist(outDir, 'dir')
    mkdir(outDir)
end
for cIdx = 1:length(maybeNRS)
    
    ftoFind = [maybeNRS(cIdx).brainArea '_' num2str(maybeNRS(cIdx).ChannelNumber) '_' num2str(maybeNRS(cIdx).Name) '_1.png'];
     
    fileNum(cIdx) = structfind(ims, 'name', ftoFind);
    
    ftoMove = [ims(fileNum(cIdx)).folder filesep ims(fileNum(cIdx)).name];
    copyfile(ftoMove, [outDir filesep ims(fileNum(cIdx)).name]);
    
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


