

% Script to help uncover yeild of face/axis tuned neurons in various patients
% read in responsive cells (could also do per session)
% Make structs with just ramp tuned neurons, max group face neurons, High
% FSI neurons
% check what patient nad side contributed to each
% vwadia May2023

% This is metadata for images showing where our targetting ended up
% relative to fMRI activation
%% Load in resp cells and find axis tuned ones

setDiskPaths
% load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])
load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'])

idx = zeros(length(strctCells), 1);

for cellIndex = l(strctCells)
    
    p_info = strctCells(cellIndex).pvalRamp;
    
    idx(cellIndex) = sum(p_info < 0.01);
end

% ramp tuned neurons
s_c = strctCells(idx >= 0.98*length(p_info));

%% set up orders

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
%% use those orders to find max group and FSI
method = 0;

% findind cell's max group
for cellIndex = l(strctCells)
    
    if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
        labels = catOrd{1};
        timelimits = [-0.17, 0.33];
        stimOffDur = 166.6250;
        stimDur = 166.6250;
        sortedOrder = order1;
    elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
        labels = catOrd{3};
        timelimits = [-0.17, 0.53];
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
    
    % find max group using smoothed psth
    n_stdDevs = 2.5;
    [~, max_group] = Utilities.computeResponseLatency(psths(cellIndex, :), labels, timelimits,...
        stimOffDur, stimDur, method, n_stdDevs);
    
    responses{cellIndex, 4} = max_group;
    
    % compute FSI
    windowBegin = (-timelimits(1)*1e3)+floor(responses{cellIndex, 2});
    windowLength = ceil(stimDur);
    faceIdx = 1; % pre-determined
    
    faceMat = psths{cellIndex, 1}(find(labels == faceIdx), :);
    nonFaceMat = psths{cellIndex, 1}(find(labels ~= faceIdx), :);
    
    if windowBegin+windowLength > size(faceMat, 2)
        faceResp = mean(mean(faceMat(:, windowBegin:end), 1));
        nonFaceResp = mean(mean(nonFaceMat(:, windowBegin:end), 1));
    else
        faceResp = mean(mean(faceMat(:, windowBegin:windowBegin+windowLength), 1));
        nonFaceResp = mean(mean(nonFaceMat(:, windowBegin:windowBegin+windowLength), 1));
    end
    FSI(cellIndex, 1) = (faceResp - nonFaceResp)/(faceResp + nonFaceResp);
    responses{cellIndex, 5} = (faceResp - nonFaceResp)/(faceResp + nonFaceResp);
    
end


f_c = strctCells(cell2mat(responses(:, 4)) == 1);
f_C = struct2cell(f_c')';
real_f_c = strctCells(cell2mat(responses(:, 5)) >= 0.33);
real_f_C = struct2cell(real_f_c')';

%% Parse out contributions of each session to this

% make cell arrays 
strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

s_C = struct2cell(s_c');
s_C = s_C';

PT = {};
pt_ids = {'P71CS', 'P73CS', 'P75CS', 'P76CS', 'P77CS_L', 'P77CS_R', 'P78CS', 'P79CS_L', 'P79CS_R', 'P80CS_L', 'P80CS_R', 'P81CS_L', 'P81CS_R', 'P82CS_L', 'P82CS_R', 'P84CS_L', 'P84CS_R', 'P85CS_L', 'P85CS_R'};
pt_ids = fliplr(pt_ids);

for i = 1:length(pt_ids)
    
  PT{i, 1} = pt_ids{i};
  
  if strcmp(pt_ids{i}(end) , 'L')
      
      t1 = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), strctCELL(:, 7)); t2 = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));
      PT{i, 2} = sum(t1+t2 == 2);
      
      t1_sc = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), s_C(:, 7)); t2_sc = cellfun(@(x) strcmp(x, 'LFFA'), s_C(:, 4));
      PT{i, 3} = sum(t1_sc+t2_sc == 2);
      
      t1_fc = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), f_C(:, 7)); t2_fc = cellfun(@(x) strcmp(x, 'LFFA'), f_C(:, 4));
      PT{i, 4} = sum(t1_fc+t2_fc == 2);
      
      t1_FSI = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), real_f_C(:, 7)); t2_FSI = cellfun(@(x) strcmp(x, 'LFFA'), real_f_C(:, 4));
      PT{i, 5} = sum(t1_FSI+t2_FSI == 2);
      
  elseif strcmp(pt_ids{i}(end) , 'R')
      
      t1 = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), strctCELL(:, 7)); t2 = cellfun(@(x) strcmp(x, 'RFFA'), strctCELL(:, 4));
      PT{i, 2} = sum(t1+t2 == 2); 
      
      t1_sc = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), s_C(:, 7)); t2_sc = cellfun(@(x) strcmp(x, 'RFFA'), s_C(:, 4));
      PT{i, 3} = sum(t1_sc+t2_sc == 2);
      
      t1_fc = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), f_C(:, 7)); t2_fc = cellfun(@(x) strcmp(x, 'RFFA'), f_C(:, 4));
      PT{i, 4} = sum(t1_fc+t2_fc == 2);
      
      t1_FSI = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), real_f_C(:, 7)); t2_FSI = cellfun(@(x) strcmp(x, 'RFFA'), real_f_C(:, 4));
      PT{i, 5} = sum(t1_FSI+t2_FSI == 2);
      
  else
      
      PT{i, 2} = sum(cellfun(@(x) strcmp(x, pt_ids{i}), strctCELL(:, 7)));
      PT{i, 3} = sum(cellfun(@(x) strcmp(x, pt_ids{i}), s_C(:, 7)));    
      PT{i, 4} = sum(cellfun(@(x) strcmp(x, pt_ids{i}), f_C(:, 7)));
      PT{i, 5} = sum(cellfun(@(x) strcmp(x, pt_ids{i}), real_f_C(:, 7)));    
  end
  
  
end  

% convert to table for ease of reading
tab = cell2table(PT, "VariableNames", ["Patient", "Responsive units", "Axis tuned units", "Max group faces", "High FSI"]);

    
    
