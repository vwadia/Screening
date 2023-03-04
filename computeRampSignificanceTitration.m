
% script to titrate down the numebr of stimuli and compute ramp significance
% Workflow:
%      Load in cells
%      Per cell:
%          Use all stim - compute ramp significance
%          Use 400, 300, 200, 100
%          Save p-values as row in matrix (cells x stim # used)
%     Compute percentages of cells with significant ramps for a given stim #

% repeat this also with all cells from 1593 screening.


%% set paths 

setDiskPaths

%% load in cells - objects
% options.screenType = 'Object';


% taskPath = 'Object_Screening';
% stimDir = '500Stimuli';
% imageIDs = [1:500]';
% 
% % stimDir = '1593Stimuli';
% % imageIDs = [1:1593]';
% 
% % load([diskPath filesep taskPath filesep 'AllITCells_500Stim_Scrn']); RespCells = 1;
% load([diskPath filesep taskPath filesep 'ITCells_500Stim_Scrn_SigRamp']); RespCells = 0;
% 
% % load params
% load([diskPath filesep 'Objectspace' filesep stimDir filesep 'params_Alexnet_fc6_500Stimuli']);


%% Faces only
options.screenType = 'Face';


taskPath = 'Face_Screening/P73CS/2000FaceScreen_Session_1_20210326';


load([diskPath filesep taskPath filesep 'PsthandResponses']); psths = screeningPsth;
load([diskPath filesep taskPath filesep 'strctCells'])

% order = load([diskPath filesep taskPath filesep 'realOrder_P73CS_FullSynthFaces.mat']);
% order = order.imageTexture - 10;

stimDir = '667Stimuli';
imageIDs = unique(order)';


load([diskPath filesep 'ObjectSpace' filesep 'parameters_2k_synthetic_faces'])
params = params(imageIDs, :);


index = cellfun(@isempty, responses);
strctCells(index(:, 1)) = [];
responses(index(:, 1), :) = [];
psths(index(:, 1), :) = [];

%% pick out visually responsive neurons

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';


if strcmp(stimDir, '500Stimuli') && RespCells
    imageIDs = [1:500]';
    order1 = repelem(imageIDs, 6);
    order2 = repelem(imageIDs, 4);
    
    
    basicMethod = 0;
    compResponses = cell(length(strctCells), 3);
    
    for cellIndex = l(strctCells)
        
        if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
            %         labels = catOrd{1};
            [labels, ~] = Utilities.makeObjCatLabsScreening('P75CS_ObjScreen', order1);
            timelimits = [-0.17, 0.33];
            stimOffDur = 166.6250;
            stimDur = 166.6250;
            sortedOrder = order1;
        else
            %         labels = catOrd{2};
            [labels, ~] = Utilities.makeObjCatLabsScreening('P75CS_ObjScreen', order2);
            timelimits = [-0.13, 0.53];
            stimOffDur = 133.4680;
            stimDur = 266.6250;
            sortedOrder = order2;
        end
        %      labels = catOrd;
        %      timelimits = [-0.13, 0.53];
        %      stimOffDur = 133.4680;
        %      stimDur = 266.6250;
        %      [sortedOrder, correctOrder] = sortrows(order);
        
        n_stdDevs = 2;
        [respLat, max_group] = Utilities.computeResponseLatency(psths(cellIndex, :), labels, timelimits,...
            stimOffDur, stimDur, basicMethod, n_stdDevs);
        
        endRas = size(psths{cellIndex, 1}, 2);
        if respLat ~= 0
            %         respLat = 100 + (-screeningData.timelimits(1)*1e3); % choose this manually
            windowLength = floor(stimDur);
            windowBegin = respLat;
            windowEnd = windowBegin+windowLength;
            if windowEnd > endRas
                windowEnd = endRas;
            end
            
            
            for i = l(imageIDs)
                stimRaster = psths{cellIndex, 1}(find(sortedOrder == imageIDs(i)), windowBegin:windowEnd);
                compResponses{cellIndex, 1}(i, 1) = mean(mean(stimRaster))*1e3; % note the 1000x multiplcation
            end
            
            compResponses{cellIndex, 2} = respLat - (-timelimits(1)*1e3);
            compResponses{cellIndex, 3} = strctCells(cellIndex).Name;
            if exist('max_group')
                compResponses{cellIndex, 4} = max_group;
            end
            respLat = 0; % empty it so it gets assigned again if present
        end
        %     if exist('max_group')
        %         compResponses{cellIndex, 4} = max_group;
        %     end
        index = cellfun(@isempty, compResponses);
        strctCells(index(:, 1)) = [];
        responses = compResponses; responses(index(:, 1), :) = [];
        psths(index(:, 1), :) = [];
    end
end




%% compute ramp significance


% ind_train = [length(imageIDs):-25:25];

ind_train = [667:-66:0];


sigMat = nan(length(strctCells), length(ind_train));

for cellIndex = l(strctCells)
    tic
    ctr = 1;
    for ind = ind_train
        
        n_reps = 50;
        % randomly select paramters
        options.ind_train = sortrows(randsample(imageIDs, ind, false));        
        p_dist = nan(1, n_reps);
        for dist = 1:n_reps
            [p, pRat] = Utilities.ObjectSpace.linearity_measure_STA(responses{cellIndex, 1}, params, options);
            p_dist(dist) = p;
        end
        
        p_sum = sum(0.01 < p_dist)/length(p_dist);
        sigMat(cellIndex, ctr) = 1-p_sum;
    
        ctr = ctr + 1;
    end
    toc
    disp(['Finished for cell ' num2str(cellIndex)])
%     keyboard
    
    
end

if strcmp(options.screenType, 'Object')
    sM = double(sigMat < 0.01);
elseif strcmp(options.screenType, 'Face')
    id = sigMat(:, 1) ~= 1;
    sM = double(sigMat);
    sM(id, :) = [];
    
end

%%
f = figure;
x = categorical(ind_train);
% b = bar(x, sum(sM)/length(strctCells));
b = bar(x, sum(sM)/length(sM));
% b = bar(x, sM);
b.FaceColor = [0.4940 0.1840 0.5560];
set(gca, 'xdir', 'reverse')
xlabel('Number of stimuli used to compute STA')
ylabel('Fraction of cells with significant ramps')

if strcmp(options.screenType, 'Object')
    title('Ramp significance as a function of number of stimuli')
    if RespCells
        filename = [diskPath filesep taskPath filesep 'RampSig_vs_NoOfStim_RespCells'];
    else
        filename = [diskPath filesep taskPath filesep 'RampSig_vs_NoOfStim_SigCells'];
    end
    
elseif strcmp(options.screenType, 'Face')
    title('Ramp significance - Face cells')
    filename = [diskPath filesep taskPath filesep 'RampSig_vs_NoOfStim_FaceCellsSAScreen'];
    
end
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% 

print(f, filename, '-dpng', '-r0')

%% Break Point analysis - point where it is not guaranteed that a cell will show up ramp tuned 100% of the time

bP = sM ~= 1;
brk_pt = nan(size(bP, 1), 1);

for i = 1:size(bP, 1)
    
    
    row = bP(i, :);
    brk = find(row == 1); brk_pt(i) = brk(1);
    
    
end
brk_pt = ind_train(brk_pt);

%% break point - histogram/heatmap

% f2 = figure; 
% histogram(brk_pt, 'FaceColor', [0.6350 0.0780 0.1840], 'EdgeColor', [0 0 0]); 
% % histogram(brk_pt,  0:66:660, 'FaceColor', [0.6350 0.0780 0.1840], 'EdgeColor', [0 0 0]); 
% filename = [diskPath filesep taskPath filesep 'BreakPoint_Ramptuning_SAScreen'];
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/length(sM), 2))
% title('Drop point analysis - Face cells - SA Screen')
% xlabel('Number of stimuli used to compute STA')
% ylabel('Fraction of cells with significant ramps')
% set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% 
% print(f2, filename, '-dpng', '-r0')
% 
% 
% close all

% heatmap
f2 = figure; 
imagesc(sM);

% set up colormap
orangemap = esa(300);
[WhiteColor, Whitepos] = max(sum(orangemap,2));
orangemap = orangemap([1:Whitepos round(linspace(Whitepos+1,size(orangemap,1)-2,Whitepos-1))],:);
colormap(orangemap)


% set up colorbar -------------------------------------------------
cb = colorbar;
cb.Ticks = ([0 1]);
cb.FontWeight ='bold';
cb.Position = [0.920833333333333,0.719047619047619,0.025595238095239,0.204761904761913];


filename = [diskPath filesep taskPath filesep 'BreakPointHeatMap_Ramptuning_SAScreen'];
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel', [667:-66:7])
title('Probability of ramp significance - Face cells')
xlabel('Number of stimuli used to compute STA')
ylabel('Cell number')



set(gca, 'FontSize', 14, 'FontWeight', 'bold')
print(f2, filename, '-dpng', '-r0')


close all



