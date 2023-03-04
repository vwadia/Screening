

setDiskPaths

taskPath = 'Object_Screening';
load([diskPath filesep 'Object_Screening' filesep 'ITCells_500stim_Scrn_SigRamp']);


%% set up orders

imageIDs = [1:500]';
order1 = repelem(imageIDs, 6);
order2 = repelem(imageIDs, 4);
catOrd = {};
% make cat labels
anovaType = 'CategoryObject';
faceInds = 134:210;
objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
textInds = [264:282 290 400:408];
vegInds = [211:235 357:399];
animInds = [1:84 256:263];

for i = 1:2

    if i == 1
        order = order1; % 3000 elements
    else
        order = order2; % 2000 elements
    end

    catOrder = zeros(length(order), 1);
    catOrder(ismember(order, faceInds)) = 1;
    catOrder(ismember(order, textInds)) = 2;
    catOrder(ismember(order, vegInds)) = 3;
    catOrder(ismember(order, animInds)) = 4;
    catOrder(ismember(order, objInds)) = 5;

    catOrd{i} = catOrder;
end


%% compute new responses


basicMethod = 0;
compResponses = cell(length(strctCells), 3);

for cellIndex = 1:length(strctCells)

    if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
        labels = catOrd{1};
        timelimits = [-0.17, 0.33];
        stimOffDur = 166.6250;
        stimDur = 166.6250;
        sortedOrder = order1;
    else
        labels = catOrd{2};
        timelimits = [-0.13, 0.53];
        stimOffDur = 133.4680;
        stimDur = 266.6250;
        sortedOrder = order2;
    end



    [respLat, max_group] = Utilities.computeResponseLatency(psths(cellIndex, :), labels, timelimits,...
        stimOffDur, stimDur, basicMethod);

    endRas = size(psths{cellIndex, 1}, 2);
    if respLat ~= 0
        %         respLat = 100 + (-screeningData.timelimits(1)*1e3); % choose this manually
        windowLength = floor(stimDur);
        windowBegin = respLat;
        windowEnd = windowBegin+windowLength;
        if windowEnd > endRas
            windowEnd = endRas;
        end


        for i = 1:length(imageIDs)
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

    if basicMethod == 0
        responses = compResponses;
    end

    strctCells_cR = strctCells;
    responses_cR = responses;
    psths_cR = psths;

end
%% save the categories and/or plot them
for cat = 1:5

    ids = find(cell2mat(compResponses(:, 4)) == cat);
    strctCells = strctCells_cR(ids);
    responses = responses_cR(ids, :);
    psths = psths_cR(ids, :);

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


    save([diskPath filesep 'Object_Screening' filesep 'splitByCategory' filesep ['IT' group '_SigRamp_500Stimuli.mat']], 'strctCells', 'responses', 'psths')
end