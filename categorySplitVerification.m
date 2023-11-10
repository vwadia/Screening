
% script to produce STA plots and psths to ensure I have divvyed up cells
% by category correctly

% load in specific category mat file
% produce psths, and STA plots




setDiskPaths

taskPath = 'Object_Screening';
% load([diskPath filesep 'Object_Screening' filesep 'ITCells_500stim_Scrn_SigRamp']);
% 
% outPath = [diskPath filesep taskPath filesep 'SigRamp'];



%% load in and plot the categories

for cat = 1:5

    switch cat
        case 1
            group = 'faceCells';
        case 2
            group = 'textCells';
        case 3
            group = 'vegCells';
        case 4
            group = 'animCells';
        case 5
            group = 'objCells';
    end
    outPath = [diskPath filesep taskPath filesep 'SigRamp' filesep group];

    load([diskPath filesep 'Object_Screening' filesep 'splitByCategory' filesep ['IT' group '_SigRamp_500Stimuli.mat']]);
% end
    subID = 'AllCells_ParamObj';

    strctCELL = struct2cell(strctCells);
    strctCELL = reshape(strctCELL, [size(strctCELL, 1) size(strctCELL, 3)]);
    strctCELL = strctCELL';
    P71Cells = cellfun(@(x) strcmp(x, 'P71CS'), strctCELL(:, 7));

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


    % plot to confirm - P71 session with different settings

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



    %%
    rest_cells = strctCells(~P71Cells); % not counting the ones from P71
    sortedOrder = repelem([1:500], 4)';
    stimDur = 266.6250;
    timelimits = [-0.13 0.53];

    [labels, anovaType] = Utilities.makeObjCatLabsScreening(subID, sortedOrder);

    % dataParams.sigCellAnnoying = dataParams.sigCell;

    dataParams.responses = responses(~P71Cells, :);
    dataParams.psth = psths(~P71Cells, :);

    dataParams = Utilities.runAnova1(labels, stimDur, timelimits, rest_cells, dataParams);

    %  plot to confirm - rest of the cells

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

end


images.dicom.parseDICOMDIR
