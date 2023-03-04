% Script to load in and compare cells from morning and afternoon sessions
% to see if they are the same
% 
% vwadia/March 2022

% Edits May2022
% 
% Cells in the morning and afternoon *must* have the same ramp status. 
% Options: Check only cells with significant ramps on a given channel/area
% Compute pvalue of ramp

%% paths to stuff, set parameters, and load in data
dbstop if error
setDiskPaths

ITCellsOnly = 1;

imageIDs = [1:500]';
stimDur = 267; 
channels = [209:224]; % differs by session
morn = struct;  aft = struct;
order = repelem([1:500]', 4);

% m_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_1_20210917']; m_subID = 'P76CSFast_2';
% a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917']; a_subID = 'P76CSRec_ReScreen';
% load([a_basePath filesep 'PsthandResponses']); aft.order = order;  order = repelem([1:500]', 4); % because order is different here - only needed for full structure way of loading data. reset order after assignment to structure

% m_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_3_20210927']; m_subID = 'P76CS_RecScreen_3';
% a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927']; a_subID = 'P76CSRec_ReScreen_3';

% m_basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'RecallScreening_Session_1_20220330']; m_subID = 'P79CS_1';
% a_basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330']; a_subID = 'P79CS_ReScreen_1';

% m_basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'RecallScreening_Session_2_20220403']; m_subID = 'P79CS_3';
% a_basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403']; a_subID = 'P79CS_ReScreen_3';

% m_basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'RecallScreening_Session_3_20220405']; m_subID = 'P79CS_4';
% a_basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405']; a_subID = 'P79CS_ReScreen_4';

% m_basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'RecallScreening_Session_1_20220728']; m_subID = 'P80CS_RecScreen_1';
% a_basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728']; a_subID = 'P80CS_ReScreenRecall';

% m_basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'RecallScreening_Session_2_20220730']; m_subID = 'P80CS_RecScreen_2';
% a_basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']; a_subID = 'P80CS_ReScreecRecall_2';

m_basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030']; m_subID = 'P81CS_AM';
a_basePath = [diskPath filesep 'Object_screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030']; a_subID = 'P81_synth';



%% load in data - all responsive cells from both sets of sessions. Will have to remake large struct for other sessions

load([diskPath filesep 'Object_screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])
withRampPValDist = 1;

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

% carve out morning session
m_c = cellfun(@(x) strcmp(x, m_subID), strctCELL(:, 8));

morn.psths = psths(m_c, :);
morn.responses = responses(m_c, :);
morn.strctCells = strctCells(m_c);
morn.order = order;

% carve out afternoon session
a_c = cellfun(@(x) strcmp(x, a_subID), strctCELL(:, 8));

aft.psths = psths(a_c, :);
aft.responses = responses(a_c, :);
aft.strctCells = strctCells(a_c);
if ~isfield(aft, 'order')
    aft.order = order;
end

SSIM = 0;

%% load in data - with full cell structures. Not as flexible for other areas

% % all cells - morning
% load([diskPath filesep 'Object_Screening' filesep 'ALLITCells_500Stim_Scrn.mat']);
% withRampPValDist = 0;

% 
% m_strctCELL = struct2cell(strctCells');
% m_strctCELL = m_strctCELL';
% 
% m_c = cellfun(@(x) strcmp(x, m_subID), m_strctCELL(:, 8));
% 
% morn.psths = psths(m_c, :);
% morn.responses = responses(m_c, :);
% morn.strctCells = strctCells(m_c);
% morn.order = order;
% 
% % get rid of non-responsive units 
% index = cellfun(@isempty, morn.responses);
% morn.responses(index(:, 1), :) = [];
% morn.psths(index(:, 1), :) = [];
% morn.strctCells(index(:, 1)) = [];
% 
% % all cells - afternoon. Note that this will reset 'responses'/'psths'/'strctCells'
% load([diskPath filesep 'Object_Screening' filesep 'ALLITCells_500Stim_ReScreen.mat']);
% 
% a_strctCELL = struct2cell(strctCells');
% a_strctCELL = a_strctCELL';
% 
% a_c = cellfun(@(x) strcmp(x, a_subID), a_strctCELL(:, 8));
% 
% aft.psths = psths(a_c, :);
% aft.responses = responses(a_c, :);
% aft.strctCells = strctCells(a_c);
% if ~isfield(aft, 'order')
%     aft.order = order;
% end
% 
% % get rid of non-responsive units 
% index = cellfun(@isempty, aft.responses);
% aft.responses(index(:, 1), :) = [];
% aft.psths(index(:, 1), :) = [];
% aft.strctCells(index(:, 1)) = [];
% 
% SSIM = 0;

%% load in data - old way individual sessions 
% withRampPValDist = 0;
% 
% % Morning -----------------------------------------------------------------
% load([m_basePath filesep 'PsthandResponses']);
% load([m_basePath filesep 'strctCells']);
% 
% strctCELL = struct2cell(strctCells');
% strctCELL = strctCELL';
% 
% if ITCellsOnly
%     IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
%     morn.strctCells = strctCells(IT_Cells);
%     morn.psths = screeningPsth(IT_Cells, :);
%     morn.responses = responses(IT_Cells, :);
%     morn.order = order;
% end
% 
% % get rid of non-responsive units for both sessions
% index = cellfun(@isempty, morn.responses);
% morn.responses(index(:, 1), :) = [];
% morn.psths(index(:, 1), :) = [];
% morn.strctCells(index(:, 1)) = [];
% 
% 
% % Afternoon ---------------------------------------------------------------
% load([a_basePath filesep 'PsthandResponses']);
% load([a_basePath filesep 'strctCells']);
% 
% strctCELL = struct2cell(strctCells');
% strctCELL = strctCELL';
% 
% if ITCellsOnly
%     IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
%     aft.strctCells = strctCells(IT_Cells);
%     aft.psths = screeningPsth(IT_Cells, :);
%     aft.responses = responses(IT_Cells, :);
%     aft.order = order;
% end
% 
% % get rid of non-responsive units for both sessions
% index = cellfun(@isempty, aft.responses);
% aft.responses(index(:, 1), :) = [];
% aft.psths(index(:, 1), :) = [];
% aft.strctCells(index(:, 1)) = [];
% 
% SSIM = 0;



%% cellMat

% correctly set timelimits morning
if size(morn.psths{1, 1}, 2) < 700
    morn.timelimits = [-0.13 0.53];
elseif size(morn.psths{1, 1}, 2) >= 700
    morn.timelimits = [-0.17 0.53];
end

% correctly set timelimits afternoon
if size(aft.psths{1, 1}, 2) < 700
    aft.timelimits = [-0.13 0.53];
elseif size(aft.psths{1, 1}, 2) >= 700
    aft.timelimits = [-0.17 0.53];
end


if ~SSIM
    cellMat = zeros(size(morn.strctCells, 2) + size(aft.strctCells, 2), length(imageIDs));
    
    % using rank order of stim
    for m_row = 1:size(morn.strctCells, 2)
        windowBegin = morn.responses{m_row, 2}-morn.timelimits(1)*1e3;
        if ~isempty(windowBegin)
            cellMat(m_row, :) = Utilities.sortByRespMagnitude(morn.order, imageIDs, morn.psths{m_row, 1}, windowBegin, stimDur)';
        else
            disp('hi')
            cellMat(m_row, :) = Utilities.sortByRespMagnitude(morn.order, imageIDs, morn.psths{m_row, 1}, [], stimDur, morn.timelimits)';           
        end
    end
    
    for a_row = 1:size(aft.strctCells, 2)
        windowBegin = aft.responses{a_row, 2}-aft.timelimits(1)*1e3;
        if ~isempty(windowBegin)
            cellMat(a_row+size(morn.strctCells, 2), :) = Utilities.sortByRespMagnitude(aft.order, imageIDs, aft.psths{a_row, 1}, windowBegin, stimDur)';
        else
            disp('hi')
            cellMat(a_row+size(morn.strctCells, 2), :) = Utilities.sortByRespMagnitude(aft.order, imageIDs, aft.psths{a_row, 1}, [], stimDur, aft.timelimits)';   
        end
    end
    
% cellMat - for SSIM
else
    sc_Psth = [morn.psths(:, 1); aft.psths(:, 1)];
    
    delIdx = [];
    for i = 1:length(imageIDs)
        if length(find(aft.order == i)) < 4
            delIdx = [delIdx; i];
        end
    end
    
    % fixing order
    numReps = 4;
    testOrd = morn.order;
    testOrd(delIdx*numReps) = [];
    assert(isequal(testOrd, aft.order), 'The order is wrong!');
    
    % adjusting morning psths
    for i = 1:size(morn.strctCells, 2)
        
        sc_Psth{i, 1}(delIdx*numReps, :) = [];
        
    end
    
    cellMat = sc_Psth;
end


%% cellID array

offset = 100; % matches the response latency function
% putting together cellID array with morning sess labeled as 1 and
% afternoon as 2
cellIDs = [cell2mat(morn.responses(:, 3)) ones(size(morn.strctCells, 2), 1)];
cellIDs = [cellIDs; [cell2mat(aft.responses(:, 3)) ones(size(aft.strctCells, 2), 1)+1]];

% response latencies
cellIDs = [cellIDs [cell2mat(morn.responses(:, 2)); cell2mat(aft.responses(:, 2))]];

FRs = [];
for i = 1:size(morn.strctCells, 2)
    
%     cellIDs(i, 3) = morn.strctCells(i).ChannelNumber;
    FRs = [FRs; mean(mean(morn.psths{i, 1}(:, -morn.timelimits*1e3:-morn.timelimits*1e3+offset)))]; % the baseline is the first bit after stimON
    
end

for j = 1:size(aft.strctCells, 2)
    
%     cellIDs(size(morn.strctCells, 2)+j, 3) = aft.strctCells(j).ChannelNumber;
    FRs = [FRs; mean(mean(aft.psths{j, 1}(:, -aft.timelimits*1e3:-aft.timelimits*1e3+offset)))];
    
end

cellIDs(:, 4) = FRs;

%% load in waveforms and compute burst indices

% read in A*_sorted_new file
% average decorrelated traces per cell
% feed them into cellIDs
waveforms = [];
waveNames = [];
burstIdx = [];
sessions = 2;
for sess = 1:sessions
    if sess == 1
        chanPath = [m_basePath filesep 'sort' filesep 'final'];
    elseif sess == 2
        chanPath = [a_basePath filesep 'sort' filesep 'final'];
    end
    cellNames = cellIDs(:, 1);
    cellNames = cellNames(find(cellIDs(:, 2) == sess)); % to prevent cells from other session creeping in
    for chan = channels

        %     chanPath = [m_basePath filesep 'sort' filesep 'final'];
        spec_chan = [chanPath filesep 'A' num2str(chan) '_sorted_new.mat'];
        if exist(spec_chan, 'file')
            load([chanPath filesep 'A' num2str(chan) '_sorted_new.mat']); % creates assignedNegative with the names and allSpikesCorrFree/newSpikesNegative with spikes
            
            % average and extract waveforms
            for cN = 1:length(cellNames)
                

                if sum(ismember(assignedNegative, cellNames(cN))) > 0 && sum(ismember(useNegative, cellNames(cN))) > 0 % to prevent non-resp cells from other channels creeping in
                    waveforms = [waveforms; mean(allSpikesCorrFree(ismember(assignedNegative, cellNames(cN)), :), 1)];
                    waveNames = [waveNames; cellNames(cN)];
                    
                    % compute burst indices
                    ISI = diff(newTimestampsNegative(ismember(assignedNegative, cellNames(cN))));
                    ISI_sec = ISI*1e-6; % convert to seconds
                    burstIdx = [burstIdx; Utilities.calcBurstIndex(ISI_sec)];
                end
            end
            
            
        end
        
    end
end

cellIDs(:, 5) = burstIdx;

assert(isequal(waveNames, cellIDs(:, 1)), 'Check cells');

%% markig left and right sides - unnecessary if you restrict search to specific channels
 
sessions = 2; % morning and afternoon
for sess = 1:sessions
    
    left = 0;
    right = 1;
    
    if sess == 1
        sess_cells = morn.strctCells;
    elseif sess == 2
        sess_cells = aft.strctCells;
    end
    
    strctCELL = struct2cell(sess_cells');
    strctCELL = strctCELL';
    
    LeftIT = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));    
    RightIT = cellfun(@(x) strcmp(x, 'RFFA'), strctCELL(:, 4));
        
    sideLabels{sess} = RightIT; % 0 means left and 1 means right
   
    
end

cellIDs(:, 6) = [sideLabels{1}; sideLabels{2}];

%% marking which channels cells came from - in case this is needed

for i = 1:size(morn.strctCells, 2)
    
    cellIDs(i, 7) = morn.strctCells(i).ChannelNumber;
    
end

for j = 1:size(aft.strctCells, 2)
    
    cellIDs(size(morn.strctCells, 2)+j, 7) = aft.strctCells(j).ChannelNumber;
    
end

%% computing the ramp pvalues of all cells and storing - only for IT neurons
if ITCellsOnly 
    if ~withRampPValDist
        load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50
        
        options.screenType = 'Object';
        options.ind_train = imageIDs; % use all objects to calculate STA
        
        
        for cellIndex = 1:length(morn.strctCells)
            tic
            n_reps = 100;
            p_dist = nan(1, n_reps);
            for dist = 1:n_reps
                [p, pRat] = Utilities.ObjectSpace.linearity_measure_STA(morn.responses{cellIndex, 1}, params, options);
                p_dist(dist) = p;
            end
            p_sum = sum(0.01 < p_dist)/length(p_dist);
            disp(mean(p_dist));
            if p_sum > 0.1
                cellIDs(cellIndex, 8) = 0.1; % not significnat
            else
                cellIDs(cellIndex, 8) = 0.0099;
            end
            disp(['Finished for morning cell ' num2str(cellIndex)])
            toc
        end
        
        for cellIndex = 1:length(aft.strctCells)
            tic
            n_reps = 100;
            p_dist = nan(1, n_reps);
            for dist = 1:n_reps
                [p, pRat] = Utilities.ObjectSpace.linearity_measure_STA(aft.responses{cellIndex, 1}, params, options);
                p_dist(dist) = p;
            end
            p_sum = sum(0.01 < p_dist)/length(p_dist);
            disp(mean(p_dist));
            if p_sum == 0
                cellIDs(cellIndex+size(morn.strctCells, 2), 8) = 0; % significant
            else         
                cellIDs(cellIndex+size(morn.strctCells, 2), 8) = 0.1; % not significnat
            end
            disp(['Finished for afternoon cell ' num2str(cellIndex)])
            toc
        end
    elseif withRampPValDist
        for cellIndex = 1:length(morn.strctCells)
            p_idx = morn.strctCells(cellIndex).pvalRamp;
            
            if sum(p_idx < 0.01) == length(p_idx)
                cellIDs(cellIndex, 8) = 0;
            else
                cellIDs(cellIndex, 8) = 0.1;
            end
            
        end
        for cellIndex = 1:length(aft.strctCells)
            p_idx = aft.strctCells(cellIndex).pvalRamp;
            if sum(p_idx < 0.01) == length(p_idx)
                cellIDs(cellIndex+size(morn.strctCells, 2), 8) = 0;
            else
                cellIDs(cellIndex+size(morn.strctCells, 2), 8) = 0.1;
            end
        end

    end
end
% keyboard
%% call compareCells function 
%% after restriction to channel - with waveform does better than without 

% separate by channel 
restrictToChannels = 1;
no_repeats = 1;

% without waveform - worse performance after channel restriction
% [cellPairs, compMat] = Utilities.compareCells(cellMat, cellIDs, restrictToChannels);

% with waveform - works better April/2022
[cellPairs, compMat] = Utilities.compareCells_2(cellMat, cellIDs, restrictToChannels, waveforms, no_repeats);

% Maybe include angle between axes? or spread of stimuli in PC1-PC2 space?
% HOw does this function even work?

%% adjustment for category - makes performance worse

% faceInds = 134:210;
% objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
% textInds = [264:282 290 400:408];
% vegInds = [211:235 357:399];
% animInds = [1:84 256:263];
% 
% t = ismember(cellMat, faceInds).*1e3 + ismember(cellMat, objInds).*2e3 + ismember(cellMat, textInds).*3e3 + ismember(cellMat, vegInds).*4e3 + ismember(cellMat, animInds).*5e3; 
% 
% cellMat = cellMat + t;

% % cellMat = zeros(length(order), 1);
% cellMat(ismember(cellMat, faceInds)) = 1;
% cellMat(ismember(cellMat, textInds)) = 2;
% cellMat(ismember(cellMat, vegInds)) = 3;
% cellMat(ismember(cellMat, animInds)) = 4;
% cellMat(ismember(cellMat, objInds)) = 5;


