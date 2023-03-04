
% script to produce STA plots and psths for all IT cells
% appropriately correcting for the annoying session of P71

% load in specific category mat file
% produce psths, and STA plots




setDiskPaths

taskPath = 'Object_Screening';
% load([diskPath filesep 'Object_Screening' filesep 'ITCells_500stim_Scrn_SigRamp']); 
% load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500stim_Scrn_SigRamp']); 


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

outPath = [diskPath filesep taskPath filesep 'forPaper' filesep 'AllCells'];
if ~exist(outPath)
    mkdir(outPath)
end

%%
% [-0.17 0.33], 166.6250, numReps = 6
% [-0.13 0.53], 266.6250, numReps = 4
% load([taskPath filesep 'AllCells_500stim_Scrn']);

% get rid of non-responsive units
% index = cellfun(@isempty, responses);
% responses(index(:, 1), :) = [];
% psths(index(:, 1), :) = [];
% strctCells(index(:, 1)) = [];


subID = 'AllCells_ParamObj';

strctCELL = struct2cell(strctCells);
strctCELL = reshape(strctCELL, [size(strctCELL, 1) size(strctCELL, 3)]);
strctCELL = strctCELL';
P71Cells = cellfun(@(x) strcmp(x, 'P71CS'), strctCELL(:, 7));
P76_ShortSessCells = cellfun(@(x) strcmp(x, 'P76CSRec_ReScreen'), strctCELL(:, 8));


annoyingSess_cells = strctCells(P71Cells); 
sortedOrder = repelem([1:500], 6)';
stimDur = 166.6250;
stimOffDur = 166.6250;
timelimits = [-0.17 0.33];

% get labels - to make category Psths. 
[labels, anovaType] = Utilities.makeObjCatLabsScreening(subID, sortedOrder);



dataParams = struct;
% responses, timelimits, psths, stimDur
dataParams.responses = responses(P71Cells, :);
dataParams.psth = psths(P71Cells, :); 

dataParams = Utilities.runAnova1(labels, stimDur, timelimits, annoyingSess_cells, dataParams);


%% plot to confirm - P71 session with different settings

dataParams.lgnd = {'Faces', 'Text', 'Plants/fruits', 'Animals', 'Objects'};
dataParams.catOrder = labels;
dataParams.catIDs = unique(labels);
dataParams.timelimits = timelimits;
dataParams.anovaType = anovaType;
dataParams.Binsize = floor(stimDur*0.1);
dataParams.stimDur = stimDur;
dataParams.stimOffDur = stimOffDur;
dataParams.subID = 'Placeholder';
imagesPerSubplot = length(dataParams.catIDs);
subPlotNum = 1;
numFigsPerCell = 1;

Utilities.Plotting.PlotRastersScreening(outPath, imagesPerSubplot, subPlotNum, numFigsPerCell, annoyingSess_cells, dataParams, 0)
disp('Done')


%% STA plots to confirm
% load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50
% pathOut = [taskPath filesep 'STA_and_projections_varWindow'];
% 
% if ~exist(pathOut)
%     mkdir([pathOut]);
% end
% dataParams.imageIDs = [1:500]';
% options.screenType = 'Object';
% 
% alpha = 0.01;
% % for cellIndex = l(annoyingSess_cells)
% for cellIndex = l(strctCells)
%     
%     if ~isempty(dataParams.responses{cellIndex, 2})
%         options.ind_train = dataParams.imageIDs; % use all objects to calculate STA
%         [hfig, p, options] = Utilities.ObjectSpace.STA_figure_original(dataParams.responses{cellIndex}, params, options); % pass score to this instead of projectedResponses
%         if isfield(options, 'xvals') && isfield(options, 'yvals')
%             strctCells(cellIndex).Im_xvals = options.xvals;
%             strctCells(cellIndex).Im_yvals = options.yvals;
%         end
%         sgtitle({['Cell number ' num2str(annoyingSess_cells(cellIndex).Name)] 'STA and projections ', annoyingSess_cells(cellIndex).brainArea});
%         annoyingSess_cells(cellIndex).pvalRamp = p;
%         
%         if p < alpha
%             newPathOut = [pathOut filesep 'significant_cells'];
%             if ~exist(newPathOut)
%                 mkdir([newPathOut]);
%             end
%             if isfield(options, 'encoded_stim') || isfield(options, 'recalled_stim')
%                 print(hfig, [newPathOut filesep annoyingSess_cells(cellIndex).brainArea '_' num2str(annoyingSess_cells(cellIndex).ChannelNumber) '_' num2str(annoyingSess_cells(cellIndex).Name) '_orderedStim'], '-dpng', '-r0')
%             else
%                 print(hfig, [newPathOut filesep annoyingSess_cells(cellIndex).brainArea '_' num2str(annoyingSess_cells(cellIndex).ChannelNumber) '_' num2str(annoyingSess_cells(cellIndex).Name)], '-dpng', '-r0')
%             end
%         else
%             print(hfig, [pathOut filesep annoyingSess_cells(cellIndex).brainArea '_' num2str(annoyingSess_cells(cellIndex).ChannelNumber) '_' num2str(annoyingSess_cells(cellIndex).Name)], '-dpng', '-r0')
%         end
% 
%      
%         close all;
%     end
% end

%% short session cells 

shortSess_cells = strctCells(P76_ShortSessCells);

% load special order
a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];
a_strct = load([a_basePath filesep 'PsthandResponses']);
sortedOrder = a_strct.order;
clearvars a_strct

stimDur = 266.6250;
timelimits = [-0.17 0.53];

[labels, anovaType] = Utilities.makeObjCatLabsScreening(subID, sortedOrder);

dataParams.responses = responses(P76_ShortSessCells, :);
dataParams.psth = psths(P76_ShortSessCells, :); 

dataParams = Utilities.runAnova1(labels, stimDur, timelimits, shortSess_cells, dataParams);


%% plot short session cells

dataParams.lgnd = {'Faces', 'Text', 'Plants/fruits', 'Animals', 'Objects'};
dataParams.catOrder = labels;
dataParams.catIDs = unique(labels);
dataParams.timelimits = timelimits;
dataParams.anovaType = anovaType;
dataParams.Binsize = floor(stimDur*0.1);
dataParams.stimDur = stimDur;
dataParams.subID = 'Placeholder';
imagesPerSubplot = length(dataParams.catIDs);
subPlotNum = 1;
numFigsPerCell = 1;

Utilities.Plotting.PlotRastersScreening(outPath, imagesPerSubplot, subPlotNum, numFigsPerCell, shortSess_cells, dataParams, 0)
disp('Done - ShortSess Cells')


%% rest of cells

rest_cells = strctCells(~(P71Cells+P76_ShortSessCells)); % not counting the ones from P71
sortedOrder = repelem([1:500], 4)';
stimDur = 266.6250;
timelimits = [-0.17 0.53];

[labels, anovaType] = Utilities.makeObjCatLabsScreening(subID, sortedOrder);

% dataParams.sigCellAnnoying = dataParams.sigCell;

dataParams.responses = responses(~(P71Cells+P76_ShortSessCells), :);
dataParams.psth = psths(~(P71Cells+P76_ShortSessCells), :); 

dataParams = Utilities.runAnova1(labels, stimDur, timelimits, rest_cells, dataParams);

%%  plot to confirm - rest of the cells
tic
dataParams.lgnd = {'Faces', 'Text', 'Plants/fruits', 'Animals', 'Objects'};
dataParams.catOrder = labels;
dataParams.catIDs = unique(labels);
dataParams.timelimits = timelimits;
dataParams.anovaType = anovaType;
dataParams.Binsize = floor(stimDur*0.1);
dataParams.stimDur = stimDur;
dataParams.subID = 'Placeholder';
imagesPerSubplot = length(dataParams.catIDs);
subPlotNum = 1;
numFigsPerCell = 1;

Utilities.Plotting.PlotRastersScreening(outPath, imagesPerSubplot, subPlotNum, numFigsPerCell, rest_cells, dataParams, 0)
disp('Done - Reg cells')
toc
%% STA plots - rest cells

% alpha = 0.01;
% for cellIndex = l(rest_cells)
%     
%     if ~isempty(dataParams.responses{cellIndex, 2})
%         options.ind_train = dataParams.imageIDs; % use all objects to calculate STA
%         [hfig, p, options] = Utilities.ObjectSpace.STA_figure_original(dataParams.responses{cellIndex}, params, options); % pass score to this instead of projectedResponses
%         if isfield(options, 'xvals') && isfield(options, 'yvals')
%             strctCells(cellIndex).Im_xvals = options.xvals;
%             strctCells(cellIndex).Im_yvals = options.yvals;
%         end
%         sgtitle({['Cell number ' num2str(rest_cells(cellIndex).Name)] 'STA and projections ', rest_cells(cellIndex).brainArea});
%         rest_cells(cellIndex).pvalRamp = p;
%         if p < alpha
%             newPathOut = [pathOut filesep 'significant_cells'];
%             if ~exist(newPathOut)
%                 mkdir([newPathOut]);
%             end
%             if isfield(options, 'encoded_stim') || isfield(options, 'recalled_stim')
%                 print(hfig, [newPathOut filesep rest_cells(cellIndex).brainArea '_' num2str(rest_cells(cellIndex).ChannelNumber) '_' num2str(rest_cells(cellIndex).Name) '_orderedStim'], '-dpng', '-r0')
%             else
%                 print(hfig, [newPathOut filesep rest_cells(cellIndex).brainArea '_' num2str(rest_cells(cellIndex).ChannelNumber) '_' num2str(rest_cells(cellIndex).Name)], '-dpng', '-r0')
%             end
%         else
%             print(hfig, [pathOut filesep rest_cells(cellIndex).brainArea '_' num2str(rest_cells(cellIndex).ChannelNumber) '_' num2str(rest_cells(cellIndex).Name)], '-dpng', '-r0')
%         end
% 
%      
%         close all;
%     end
% end

