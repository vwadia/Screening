
% Script to compute depth of selectivity (DOS)
% for all neurons (MTL and IT)

% DOS = n - (sum(r_i) / r_max) / (n - 1)
% where r_i is the vector of mean responses to all categories and r_max is the max mean response or max(r_i)

%% load in cells

setDiskPaths

taskPath = 'Object_Screening';
stimDir = '500Stimuli';


ITCellsOnly = 1;

if ITCellsOnly
    
    % only load responsive cells
    % load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])
    load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'])
    
    
else
    
    % no P71 cells so no need to do both orders etc.
    load([diskPath filesep taskPath filesep 'AllCells_500Stim_Scrn']) %#ok<UNRCH>
    
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    IT_MTL_Cells = cellfun(@(x) ismember(x, {'LA', 'LH', 'RA', 'RH', 'RFFA', 'LFFA'}), strctCELL(:, 4), 'UniformOutput', false);
    strctCells = strctCells(cell2mat(IT_MTL_Cells));
    psths = psths(cell2mat(IT_MTL_Cells), :);
    responses = responses(cell2mat(IT_MTL_Cells), :);
    
end

% get rid of non-responsive units
index = cellfun(@isempty, responses);
responses(index(:, 1), :) = [];
psths(index(:, 1), :) = [];
strctCells(index(:, 1)) = [];
%% set up orders 

 imageIDs = [1:500]';
order1 = repelem(imageIDs, 6);
order2 = repelem(imageIDs, 4);
order3 = load([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917' filesep 'PsthandResponses']); % had a few trials cut off
order3 = order3.order;

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
    else
        order = order3; % 1892 elements

    end
    
    catOrder = zeros(length(order), 1);
    catOrder(ismember(order, faceInds)) = 1;
    catOrder(ismember(order, textInds)) = 2;
    catOrder(ismember(order, vegInds)) = 3;
    catOrder(ismember(order, animInds)) = 4;
    catOrder(ismember(order, objInds)) = 5;
    
    catOrd{i} = catOrder;
end

%% compute DOS for each cell


DOS = nan(length(strctCells), 1);

for cellIndex = 1:length(strctCells)
    
    
     if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
        labels = catOrd{1};
        timelimits = [-0.17, 0.33];
        stimOffDur = 166.6250;
        stimDur = 166.6250;
        sortedOrder = order1;
    else
        labels = catOrd{2};
        timelimits = [-0.17, 0.53];
        stimOffDur = 133.4680;
        stimDur = 266.6250;
        sortedOrder = order2;
     end
    
     if strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
         labels = catOrd{3};  
     end
     
    N = length(unique(labels));
    respLat = responses{cellIndex, 2};
    R_i = zeros(1, N);
    
    for ct = 1:N
        
        if respLat+ceil(stimDur) <= size(psths{cellIndex, 1}, 2)
            catPsth = psths{cellIndex, 1}(labels == ct, respLat:respLat+ceil(stimDur));
        else
            catPsth = psths{cellIndex, 1}(labels == ct, respLat:end);
        end
    
        R_i(ct) = mean(mean(catPsth));
        
    end
    
    DOS(cellIndex, 1) = (N - sum(R_i)/max(R_i)) / (N-1);
    
end


%% separate areas 

if ~ITCellsOnly
    
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    Hipp_Cells = cellfun(@(x) ismember(x, {'LH', 'RH'}), strctCELL(:, 4), 'UniformOutput', false);
    Amyg_Cells = cellfun(@(x) ismember(x, {'LA', 'RA'}), strctCELL(:, 4), 'UniformOutput', false);
    IT_Cells = cellfun(@(x) ismember(x, {'LFFA', 'RFFA'}), strctCELL(:, 4), 'UniformOutput', false);
    
    
    Hipp_DOS = DOS(cell2mat(Hipp_Cells));
    Amyg_DOS = DOS(cell2mat(Amyg_Cells));
    
    MTL_DOS = DOS(logical(cell2mat(Hipp_Cells)+cell2mat(Amyg_Cells)));
    IT_DOS = DOS(cell2mat(IT_Cells));
    
else
    IT_DOS = DOS;
end

%% histogram

cols = Utilities.distinguishable_colors(3);

f = figure; 
hold on

% histogram(IT_DOS, 42, 'FaceColor', cols(1, :));
% histogram(Hipp_DOS, 42, 'FaceColor', cols(2, :));
% histogram(Amyg_DOS, 42, 'FaceColor', cols(3, :));


% histogram(IT_DOS, 10, 'FaceColor', [0.6350 0.0780 0.1840]);
histogram(IT_DOS, 10, 'FaceColor', [0.4940 0.1840 0.5560]);
% histogram(MTL_DOS, 10, 'FaceColor', [0.4940 0.1840 0.5560]);
title('DOS for visually responsive IT neurons');
ylabel('Number of Cells');
xlabel('Depth of Selectivity');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
filename = [diskPath filesep taskPath filesep 'DOS_ITneurons'];
print(f, filename, '-dpng', '-r0');


%% DOS for (most) recorded cells - 500stim screening


