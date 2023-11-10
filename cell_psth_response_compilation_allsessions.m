% Script to update the PsthandResponses and structCells
% in all session folders after new RL computation

% this is basically the relevant parts of screeningScript_main in a loop
% I should maybe re-write the screening script and make it a package...
% March 2023 vwadia

setDiskPaths

addpath(genpath([diskPath filesep 'Code' filesep 'osortTextUI']));
% addpath(genpath('ObjectSpace'));
% addpath(genpath('synthetic_face_generator'));


%% compile massive list of sessions
%% note these are divided into sessID_1 and sessID_2 simplyu because I am lazy


sessID_1 = {['Object_Screening' filesep 'P71CS' filesep 'FastObjectScreening_Session_1_20201125'],...%d
    ['Object_Screening' filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328'],...% note I'm changing this here (used to have '500Stim' appended 3/17/23 vwadia
    ['Object_Screening' filesep 'P75CS' filesep 'FingerprintScreening_Session_1_20210902'],...%d
    ['Object_Screening' filesep 'P76CS' filesep 'FingerprintScreening_Session_1_20210915'],...%d
    ['Recall_task' filesep 'P76CS' filesep 'RecallScreening_Session_1_20210917'],...%d
    ['Recall_task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925'],...%d
    ['Recall_task' filesep 'P76CS' filesep 'RecallScreening_Session_3_20210927'],...%d
    ['Object_Screening' filesep 'P77CS' filesep 'FingerprintScreening_Session_1_20211007'],...%d
    ['Object_Screening' filesep 'P78CS' filesep 'FingerprintScreening_Session_1_20220309'],...%d
    ['Object_Screening' filesep 'P78CS' filesep 'FingerprintScreening_Session_2_20220314'],...%d
    ['Recall_Task' filesep 'P79CS' filesep 'RecallScreening_Session_1_20220330'],...
    ['Object_Screening' filesep 'P79CS' filesep 'FingerprintScreening_Session_1_20220331'],...%d
    ['Recall_task' filesep 'P79CS' filesep 'RecallScreening_Session_2_20220403'],...%d
    ['Recall_task' filesep 'P79CS' filesep 'RecallScreening_Session_3_20220405'],...
    ['Object_Screening' filesep 'P80CS' filesep 'FingerprintScreening_Session_1_20220727'],...%d
    ['Recall_task' filesep 'P80CS' filesep 'RecallScreening_Session_1_20220728'],...%d
    ['Recall_task' filesep 'P80CS' filesep 'RecallScreening_Session_2_20220730'],...
    ['Object_Screening' filesep 'P80CS' filesep 'FingerprintScreening_Session_2_20220801'],...
    ['Object_Screening' filesep 'P81CS' filesep 'FingerprintScreening_Session_1_20221026'],...
    ['Object_Screening' filesep 'P81CS' filesep 'FingerprintScreening_Session_2_20221028'],...
    ['Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030'], ...
    ['Object_Screening' filesep 'P82CS' filesep 'FingerprintScreening_Session_1_20230111'],...
    ['Object_Screening' filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115']}; %d

tStrct_1 = {'P71CS_Fast_Sub_6_Block',...
    'P73CS_ParamObj_Sub_4_Block',...
    'P75CS_ObjScreen_Sub_4_Block',...
    'P76CSFast_Sub_4_Block',...
    'P76CSFast_2_Sub_4_Block',...
    'P76CS_RecScreen_2_Sub_4_Block',...
    'P76CS_RecScreen_3_Sub_4_Block',...
    'P77CS_1_Sub_4_Block',...
    'P78CS_1_Sub_4_Block',...
    'P78CS_Screen2_Sub_4_Block',...
    'P79CS_1_Sub_4_Block',...
    'P79CS_2_Sub_4_Block',...
    'P79CS_3_Sub_4_Block',...
    'P79CS_4_Sub_4_Block',...
    'P80CS_2_Sub_4_Block',...
    'P80CS_RecScreen_1_Sub_4_Block',...
    'P80CS_RecScreen_2_Sub_4_Block',...
    'P80CS_2_Att2_Sub_4_Block',...
    '81CS_forReal_Sub_4_Block',...
    'P81CS_2_Sub_4_Block',...
    'P81CS_AM_Sub_4_Block',...
    'P82CS_1_Sub_4_Block',...
    'P82CS_CL_1_Sub_4_Block'};



sessID_2 = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
    ['Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925'],... % moved ITResponses into screening folder for sess 2
    ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
    ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
    ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731'],...
    ['Object_Screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030'], ...
    ['Object_Screening' filesep 'P82CS' filesep 'ClosedLoopReScreen_Session_1_20230115']};

tStrct_2 = {'P76CSRec_ReScreen_Sub_4_Block',...
    'P76CS_RecScreen_2_Sub_4_Block',...
    'P76CSRec_ReScreen_3_Sub_4_Block',...
    'P79CS_ReScreen_1_Sub_4_Block',...
    'P79CS_ReScreen_3_Sub_4_Block',...
    'P79CS_ReScreen_4_Sub_4_Block',...
    'P80CS_ReScreenRecall_Sub_4_Block',...
    'P80CS_ReScreecRecall_2_Sub_4_Block',...
    'P81_synth_Sub_4_Block',...
    'P82CS_CLReScreen_Sub_4_Block'};



sessID = cat(2, sessID_1, sessID_2);
tStrct = cat(2, tStrct_1, tStrct_2);

sortPath = [filesep 'sort'];
finalPath = [filesep 'final'];
rawPath = [filesep 'raw'];

% grab TTL Values
pathTaskCode = [boxPath filesep 'epicScreeningVarun'];
addpath(pathTaskCode);

% CLSynthScreen = 1; noNatSort = 0;


%% save all cells (from all regions) for all sessions
siphonOff500 = 1;

for ss = 6%1:length(sessID)
    tic 
    setTTLCodes;

    basePath = [diskPath filesep sessID{ss}];
    taskStruct = load([basePath filesep tStrct{ss}]);
    % events
    events = getRawTTLs([basePath filesep rawPath filesep 'Events.nev'], 1);
    
    % rigorous way to find patient ID
    slash_pos = strfind(basePath, filesep); % all the slash positions
    p_pos = strfind(basePath, 'P');  % find all 'P's in the name
    s_pos = find(ismember(slash_pos, p_pos-1)); % single out the P that begins a folder name
    patID = basePath(slash_pos(s_pos)+1:slash_pos(s_pos+1)-1); % grab that folder name

        % extract cells
    [strctCells, dupCells] = Utilities.extractCells([basePath filesep sortPath filesep finalPath] , basePath, [], patID, taskStruct.subID);
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    pathStimuli = [basePath filesep 'stimuliUsed'];

    if ~isempty(strfind(basePath, 'ClosedLoopReScreen'))
        CLSynthScreen = 1;
        if strcmp(patID, 'P81CS')
            noNatSort = 1;
        elseif strcmp(patID, 'P82CS')
            noNatSort = 0;
        end     
    end
    
    
    % assign log file if there is one
    lgFiles = Utilities.readInFiles(basePath, 'txt');
    if isempty(lgFiles)
        clearvars log
    else
        if length(lgFiles) > 1
            for i = 1:length(lgFiles)
%                 if strcmp(lgFiles(i).name(1:length(patID)), patID)
                if ~isempty(strfind(lgFiles(i).name, 'Screening'))
                    log = [lgFiles(i).folder filesep lgFiles(i).name];
                end
            end
        else
            log = [lgFiles.folder filesep lgFiles.name];
        end
    end
    
    
    screeningData = struct;
    screeningData.subID = taskStruct.subID;
    screeningData.psth = cell(length(strctCells), 3);
    screeningData.responses = cell(length(strctCells), 1);
    
    experimentBegin = find(events(:, 2) == EXPERIMENT_ON);
    experimentEnd = find(events(:, 2) == EXPERIMENT_OFF);
    
    
    % load events
    if exist('shortLatRun') || exist('reScreenRun') % P73CS cat screen first run or AIC rescreen
        screeningData.events = events(experimentBegin(1):experimentEnd(1), :); % first run of P73 3/25
    elseif strcmp(taskStruct.subID, 'P76CSRec_ReScreen_3') || exist('longLatRun') % P73CS second run (same session)
        screeningData.events = events(experimentBegin(2):experimentEnd(2), :); % second run of P73 3/25
    elseif strcmp(taskStruct.subID, 'P76CSRec_ReScreen')
        screeningData.events = events(experimentBegin(1):experimentEnd(1), :);
    elseif strcmp(taskStruct.subID, 'P77CS_1')
        screeningData.events = events(experimentBegin(3):experimentEnd(2), :);
    elseif strcmp(taskStruct.subID, 'P79CS_ReScreen_1')
        screeningData.events = events(experimentBegin(3):experimentEnd(2), :);
    elseif strcmp(taskStruct.subID, 'P79CS_ReScreen_3')
        screeningData.events = events(experimentBegin(2):experimentEnd(2), :);
    elseif strcmp(taskStruct.subID, 'P79CS_ReScreen_4')
        screeningData.events = events(experimentBegin(3):experimentEnd(3), :);
    elseif strcmp(taskStruct.subID, 'P80CS_2')
        screeningData.events = events(experimentBegin(2):experimentEnd(2), :);
    elseif strcmp(taskStruct.subID, 'P80CS_ReScreenRecall')
        screeningData.events = events(experimentBegin(2):experimentEnd(2), :);
    elseif strcmp(taskStruct.subID, 'P80CS_ReScreecRecall_2')
        screeningData.events = events(experimentBegin(2):experimentEnd(2), :);
    else
        screeningData.events = events(experimentBegin:experimentEnd, :);
    end
    
    if strcmp(taskStruct.subID, '81CS_forReal')
        IMAGE_ON = 4; IMAGE_OFF = 5;
        IMAGE_ON_LOG = 20; IMAGE_OFF_LOG = 21;
    end
    
    %% find image on/off times
    if ~strcmp(taskStruct.subID, '62')
        imageOffPoints = find(screeningData.events(:, 2) == IMAGE_OFF);
        imageOnPoints = find(screeningData.events(:, 2) == IMAGE_ON);
    elseif strcmp(taskStruct.subID, '62')
        imageOffPoints = find(screeningData.events(:, 2) == 20);
        imageOnPoints = find(screeningData.events(:, 2) == 15);
    end
    
    imageOffTimes = screeningData.events(imageOffPoints, 1);
    imageOnTimes = screeningData.events(imageOnPoints, 1);
    
    if isequal(length(imageOffTimes), length(imageOnTimes))
        stimDur = median(imageOffTimes - imageOnTimes)*1e-3;
        stimOffDur = (imageOnTimes(2:end) - imageOffTimes(1:end-1));
        stimOffDur = median(stimOffDur(stimOffDur < 1e6))*1e-3;
    else
        stimDur = (imageOffTimes(1) - imageOnTimes(1))*1e-3;
        stimOffDur = (imageOnTimes(2) - imageOffTimes(1))*1e-3;
    end
    screeningData.stimDur = stimDur;
    screeningData.stimOffDur = stimOffDur;
    
    
    
    %% set up order
    if ~strcmp(taskStruct.subID, '62')
        order = taskStruct.order;%repmat(taskStruct.order, [10, 1]);
        order = reshape(order, [], 1); % stacking all the columns on top of each other
    else
        order = repmat([1:1593], 1, 2)';
    end
    
    
    if exist('log', 'var')
        
        % note that this gets rid of spaces in names (relevant to closed loop
        % screen with grid of images)
        fileID = fopen(log, 'r');
        text = fscanf(fileID, '%s');
        fclose(fileID);
        
        colons = strfind(text, ';');
        
        txtVals = {};
        ctr = 1;
        for i = 1:length(colons)
            if mod(i, 3) == 1
                if i == 1
                    strpos = 1:colons(i)-1;
                    txtVals{ctr, 1} = str2num(text(strpos));
                else
                    strpos = colons(i-1)+1:colons(i)-1;
                    txtVals{ctr, 1} = str2num(text(strpos));
                end
            elseif mod(i, 3) == 2
                strpos = colons(i-1)+1:colons(i)-1;
                txtVals{ctr, 2} = str2num(text(strpos));
                
            elseif mod(i, 3) == 0
                strpos = colons(i-1)+1:colons(i)-1;
                txtVals{ctr, 3} = text(strpos);
                ctr = ctr+1;
            end
        end
        
        if strcmp(taskStruct.subID, '81CS_forReal')
            logOrder = txtVals(cell2mat(txtVals(:, 2)) == IMAGE_ON_LOG, 3);
        else
            logOrder = txtVals(cell2mat(txtVals(:, 2)) == IMAGE_ON, 3);
        end
        logOrder = cellfun(@(x) x(1:end-4), logOrder, 'UniformOutput', false);

        imDir = Utilities.readInFiles(pathStimuli);
        %     imDir = dir(fullfile(pathStimuli));
        %     imDir = imDir(~ismember({imDir.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));
        [~, natIdx] = natsortfiles({imDir.name});
        names = {};
        for i = l(imDir)
            %             name = imDir(natIdx(i)).name;
            name = imDir(i).name;
            names{i, 1} = name(1:end-4);
        end
        assert(~isempty(names), 'Add stimulus images!');
        neworder = [];
        for i = l(logOrder)
            toFind = logOrder(i);
            comp = cellfun(@(x) strcmp(x, toFind), names);
            
            neworder(i, 1) = find(comp == 1);
        end
        if strcmp(taskStruct.subID, 'P76CSRec_ReScreen')
            assert(isequal(neworder, order(1:length(logOrder))));
            order = order(1:length(logOrder));            
        else
            assert(isequal(neworder, order));
        end
        
    end
 
   
    if strcmp(taskStruct.subID, 'P73CS_ParamObj')
        load([basePath filesep 'realOrder_P73CS_ParamObj.mat']);
        order = imageTexture - 10; % image textures are labelled by the picture num + 10
        order = order';
    elseif strcmp(taskStruct.subID, 'P73CS_Full')
        load([basePath filesep 'realOrder_P73CS_FullSynthFaces.mat']);
        order = imageTexture - 10; % image textures are labelled by the picture num + 10
        order = order';
    end
    
    if ~exist('CLSynthScreen', 'var') || CLSynthScreen == 0
        screeningData.imageIDs = unique(order);
    else
        screeningData.imageIDs = unique(order);
        newIms = cellfun(@(x) strcmp(x(1), 'c'), names); % so filenames have to start with a 'c'
        screeningData.imageIDs = screeningData.imageIDs(~newIms);
    end
    screeningData.numRepetitions = taskStruct.num_blocks;
    
    %% for raster extraction and plotting
    % ------------------------------------------------------------------------
    if stimOffDur < stimDur
        tl_1 = round((stimOffDur)*1e-3, 2); %'tee el dash one'
    else
        tl_1 = round((stimDur)*1e-3, 2); %'tee el dash one'
    end
    if tl_1 == 0.13
        tl_1 = 0.17;
    end
    tl_2 = round((stimDur)*2e-3, 2);
    screeningData.timelimits = [-tl_1 tl_2]; % -0.25 0.5 for IT screens, -1 2 for others
    screeningData.Binsize = floor(stimDur*0.1); % ms, for smoothing psths
    % -------------------------------------------------------------------------
    
    [sortedOrder, correctOrder] = sortrows(order);
    screeningData.sortedOrder = sortedOrder;
    screeningData.correctOrder = correctOrder;
    % using Jan's functions
    for cellIndex = 1:length(strctCells)
        
        [psth, psth1, times] = Utilities.cellfile2rasterVarun(strctCells(cellIndex),screeningData.timelimits,screeningData.Binsize,IMAGE_ON,screeningData.events);
        screeningData.psth{cellIndex, 1} = psth1(:, correctOrder)'; % rasters
        screeningData.psth{cellIndex, 2} = psth(:, correctOrder)'; % arranges the psth so all images are clumped together
        screeningData.psth{cellIndex, 3} = times; % a linspace from timecourse(1) to timecourse(2)
    end
    
    
    %% Siphoning off the 500 common stimuli from the 1593 full obj screen
     if strcmp(taskStruct.subID, 'P73CS_ParamObj') && siphonOff500
        taskStruct.subID = 'P73CS_ParamObj_500';
        basePath = [basePath filesep '500Stim'];
        if ~exist(basePath)
            mkdir(basePath);
        end
        path500Stimuli = [diskPath filesep 'Object_Screening' filesep '500Stimuli'];
        imDir = dir(fullfile(path500Stimuli));
        imDir = imDir(~ismember({imDir.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));
        [~, natIdx] = natsortfiles({imDir.name});
        
        idx = [];
        for i = l(natIdx)
            id = imDir(i).name;
            idx(end+1, 1) = str2double(id(1:4));
        end
        
        ordOverlap = ismember(sortedOrder, idx);
        sorted500Img = sortedOrder(ordOverlap);
        
        t_rem = [];
        for i = 1:length(sorted500Img)
            t = find(sorted500Img == i);
            if length(t) > screeningData.numRepetitions
                t_rem = [t_rem; t(end)]; % should I remove the last rep or some random one in the middle?
            end
        end
        
        for cellIndex = l(strctCells)
            screeningData.psth{cellIndex, 1} = screeningData.psth{cellIndex, 1}(ismember(sortedOrder, idx), :);
            screeningData.psth{cellIndex, 1}(t_rem, :) = [];
            screeningData.psth{cellIndex, 2} = screeningData.psth{cellIndex, 2}(ismember(sortedOrder, idx), :);
            screeningData.psth{cellIndex, 2}(t_rem, :) = [];
        end
        
        % so that anova labels are made correctly
        screeningData.imageIDs = [1:500]';
        screeningData.sortedOrder = repelem(screeningData.imageIDs, screeningData.numRepetitions);
        
    end
    
    if exist('CLSynthScreen', 'var') && CLSynthScreen == 1
        
        % imageIDs is only the normal screening images
        num_synth_trials = mod(length(order), length(screeningData.imageIDs)*screeningData.numRepetitions);
        
        % fix order (natsort wasn't used during the task)
        if exist('noNatSort', 'var') && noNatSort == 1
            
            taskImDir = repelem(imDir, 4); % the way the order is currently ( without natsort)
            [~, nI] = natsortfiles(taskImDir);
            
            % sort all trials
            screeningData.psth(:, 1:2) = cellfun(@(x) x(nI, :), screeningData.psth(:, 1:2), 'UniformOutput', false);
            
        end
        
        % separate the synthetic images
        screeningData.synthPsth = cellfun(@(x) x(end-num_synth_trials+1:end, :), screeningData.psth(:, 1:2), 'UniformOutput', false);
        screeningData.synthPsth(:, 3) = screeningData.psth(:, 3);
        
        % remove them from original pile
        screeningData.psth(:, 1:2) = cellfun(@(x) x(1:end-num_synth_trials, :), screeningData.psth(:, 1:2), 'UniformOutput', false);
        
    end
    %% category labeling
    [labels, anovaType] = Utilities.makeObjCatLabsScreening(taskStruct.subID, screeningData.sortedOrder);
    labels(labels == 0) = [];
    
    if exist('catOrder')
        labels = catOrder;
    end
    
    %% compute response latency via Poisson
    method = 3;
    
    if ~exist('CLSynthScreen', 'var') || CLSynthScreen == 0
        sortOrd = screeningData.sortedOrder;
    else
        % only the real images
        sortOrd = screeningData.sortedOrder(1:(length(screeningData.imageIDs)*screeningData.numRepetitions));
    end
    
    
    for cellIndex = l(strctCells)
        
        switch method
            case 0
                n_stdDevs = 2.5;
                [respLat, max_group] = Utilities.computeResponseLatency(screeningData.psth(cellIndex, :), labels, screeningData.timelimits,...
                    screeningData.stimOffDur, screeningData.stimDur, method, n_stdDevs);
            case 1
                [respLat, ~] = Utilities.computeResponseLatency(screeningData.psth(cellIndex, :), labels, screeningData.timelimits,...
                    screeningData.stimOffDur, screeningData.stimDur); %#ok<UNRCH>
            case 2
                [respLat, ~] = Utilities.computeResponseLatency(screeningData.psth(cellIndex, :), labels, screeningData.timelimits,...
                    screeningData.stimOffDur, screeningData.stimDur); %#ok<UNRCH>
            case 3
                [respLat, ~] = Utilities.computeRespLatPoisson(screeningData.psth(cellIndex, :), labels,...
                    sortOrd, screeningData.timelimits, screeningData.stimDur, true);
                respLat = -screeningData.timelimits(1)*1e3 + respLat; % adjust to be of the same for mas the other methods
        end
        
        % can get rid of this soon
        if method ~= 3
            % manually adding respLat for closed loop screen
            if strcmp(taskStruct.subID, 'P82CS_CLReScreen')
                if strctCells(cellIndex).Name == 2360
                    respLat = 300;
                end
            elseif strcmp(taskStruct.subID, 'P82CS_CL_1')
                if strctCells(cellIndex).Name == 904
                    respLat = 350;
                elseif  strctCells(cellIndex).Name  == 449
                    respLat = 360;
                end
            end
        end
        endRas = size(screeningData.psth{cellIndex, 1}, 2);
        
        if (respLat ~= 0) && ~isnan(respLat)
            %         respLat = 100 + (-screeningData.timelimits(1)*1e3); % choose this manually
            windowLength = floor(stimDur);
            windowBegin = floor(respLat);
            windowEnd = windowBegin+windowLength;
            if windowEnd > endRas
                windowEnd = endRas;
            end
            

            
            for i = l(screeningData.imageIDs)
                stimRaster = screeningData.psth{cellIndex, 1}(find(screeningData.sortedOrder == screeningData.imageIDs(i)), windowBegin:windowEnd);
                screeningData.responses{cellIndex, 1}(i, 1) = mean(mean(stimRaster))*1e3; % note the 1000x multiplcation
            end
            
            screeningData.responses{cellIndex, 2} = respLat - (-screeningData.timelimits(1)*1e3);
            screeningData.responses{cellIndex, 3} = strctCells(cellIndex).Name;
            if exist('max_group')
                screeningData.responses{cellIndex, 4} = max_group;
            end
            
            if method ~= 3
                % trimming long response latencies - if the latency is essentially in the next trial that is
                % not likely real/useful
                if stimDur < 200 || strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
                    if screeningData.responses{cellIndex, 2} > 250
                        screeningData.responses(cellIndex, :) = [];
                    end
                elseif stimDur < 300 && screeningData.responses{cellIndex, 2} > 350
                    screeningData.responses(cellIndex, :) = [];
                end
                
            end
            if exist('CLSynthScreen', 'var') && CLSynthScreen == 1
                
                % compute responses to synthetic images
                synthOrd_postSort = repelem(1:sum(newIms), screeningData.numRepetitions)';
                for i = 1:sum(newIms)
                    synthRaster = screeningData.synthPsth{cellIndex, 1}(find(synthOrd_postSort == i), windowBegin:windowEnd);
                    screeningData.synthResponses{cellIndex, 1}(i, 1) = mean(mean(synthRaster))*1e3; % note the 1000x multiplcation
                end
                
                screeningData.synthResponses{cellIndex, 2} = respLat - (-screeningData.timelimits(1)*1e3);
                screeningData.synthResponses{cellIndex, 3} = strctCells(cellIndex).Name;
            end
            
            if method ~= 3
                respLat = 0; % empty it so it gets assigned again if present
            end
        else
            screeningData.responses{cellIndex, 1} = nan;
            screeningData.responses{cellIndex, 2} = nan;
            screeningData.responses{cellIndex, 3} = strctCells(cellIndex).Name;
            
            if exist('CLSynthScreen', 'var') && CLSynthScreen == 1
                screeningData.synthResponses{cellIndex, 1} = nan;
                screeningData.synthResponses{cellIndex, 2} = nan;
                screeningData.synthResponses{cellIndex, 3} = strctCells(cellIndex).Name;
            end
        end 
        
    end
    
    
    assert(isequal(length(strctCells), length(screeningData.responses)), 'Manually fix length/order of cells')
  
    
   
    % paranoia check
    load([basePath filesep 'PsthandResponses']);
    test_cells = load([basePath filesep 'strctCells']);
    if isequal(length(test_cells.strctCells), length(strctCells))
        if exist('screeningPsth', 'var')
            assert(isequal(screeningPsth, screeningData.psth), 'Psths are different - check')
        elseif exist('psths', 'var')
            assert(isequal(psths, screeningData.psth), 'Psths are different - check')
        end
        assert(isequal(order, sortOrd), 'Orders do not match - check')
        assert(isequal(test_cells.strctCells, strctCells), 'Different cells being extracted')
    end
    clearvars psths screeningPsth responses order % remove after checks
    
    %% ANOVA
    if ~strcmp(anovaType, '1xN')
        screeningData.catOrder = labels;
        screeningData.catIDs = unique(labels);
    end
    screeningData = Utilities.selectivityCheck_ANOVA(labels, anovaType, strctCells, screeningData);
    screeningData.anovaType = anovaType;
    
    responses = screeningData.responses;
    psths = screeningData.psth;
    
    if ~exist('CLSynthScreen', 'var') || CLSynthScreen == 0
        order = screeningData.sortedOrder;
    else
        % only the real images
        order = screeningData.sortedOrder(1:(length(screeningData.imageIDs)*screeningData.numRepetitions));
    end
    
    %% save   
    if exist('IT_MTL_Cells')
        save([basePath filesep 'IT_MTL_PsthandResponses'], 'psths', 'responses', 'order', '-v7.3')
        save([basePath filesep 'IT_MTLCells'], 'strctCells')        
    else
        save([basePath filesep 'PsthandResponses'], 'psths', 'responses', 'order', '-v7.3')
        save([basePath filesep 'strctCells'], 'strctCells')
        
    end
    
    if exist('CLSynthScreen', 'var') && CLSynthScreen == 1
        synthResponses = screeningData.synthResponses;
        synthPsths = screeningData.synthPsth;
        synthOrder = synthOrd_postSort;
        
        save([basePath filesep 'SynthPsthandResponses'], 'synthPsths', 'synthResponses', 'synthOrder', '-v7.3')
    end
    
    %}
    lastsize = fprintf('Finished for session %d of %d ',ss,length(sessID));
    disp([basePath ': ' patID])
    clearvars psths responses order strctCells
toc
end





