%% script to extract responses to images during screening
%% This is a cleaned up version of the original script
%% 'screeningResponseVarun'. use this instead
%% vwadia/Jan2021


%% paths to stuff and extract cells
setDiskPaths


addpath(genpath([diskPath filesep 'Code' filesep 'osortTextUI']));
addpath(genpath('ObjectSpace'));
addpath(genpath('synthetic_face_generator'));


% basePath = [diskPath filesep 'Object_Screening' filesep 'P62CS' filesep 'RapidObjScreen_Session_1_20190417'];
% taskStruct = load([basePath filesep 'P62CS_Sub_1_Block']); patID = 'P62CS';
% FFAChansOnly = 1;

% basePath = 'Z:\dataRawEpilepsy\P69CS\10272020_screeningVarun';
% taskStruct = load([basePath filesep 'P69CS_2_Sub_6_Block']);

% basePath = 'Z:\dataRawEpilepsy\P70CS\04112020_varunScreen';
% taskStruct = load([basePath filesep 'P70CS_3_Sub_6_Block']);

% basePath = 'Z:\dataRawEpilepsy\P70CS\07112020_varunScreen';
% taskStruct = load([basePath filesep 'P70CS_4_Sub_6_Block']);

% basePath = 'Z:\dataRawEpilepsy\P70CS\11122020_varunScreen';
% taskStruct = load([basePath filesep 'P70CSDay2_Sub_6_Block']);

% basePath = [diskPath filesep 'Object_Screening' filesep 'P71CS' filesep 'CatScreen_Session_1_20201118'];
% taskStruct = load([basePath filesep 'P71CS_Sub_6_Block']); patID = 'P71CS';

% basePath = [diskPath filesep 'Recall_Task' filesep 'P71CS' filesep 'ObjectScreening_Session_1_20201121'];
% taskStruct = load([basePath filesep 'P71CS_Object_Sub_6_Block']); patID = 'P71CS'; % Nov 21
% FFAChansOnly = 1;

% basePath =  [diskPath filesep 'Object_Screening' filesep 'P71CS' filesep 'LargeObjectScreening_Session_1_20201123'];
% taskStruct = load([basePath filesep 'P71CS_LargeObject_Sub_6_Block']); patID = 'P71CS';
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P71CS' filesep 'ObjectScreening_Session_2_20201124'];
% taskStruct = load([basePath filesep 'P71CS_RecScreen2_Sub_6_Block']); patID = 'P71CS';

% basePath = [diskPath filesep 'Object_Screening' filesep 'P71CS' filesep 'FastObjectScreening_Session_1_20201125'];
% taskStruct = load([basePath filesep 'P71CS_Fast_Sub_6_Block']); patID = 'P71CS';
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P73CS' filesep 'CatScreen_Session_1_20210325'];
% % Did 2 screening runs in the same session with P73CS --------------------
% taskStruct = load([basePath filesep 'P73CS_Sub_6_Block']); shortLatRun = 1;
% % taskStruct = load([basePath filesep 'P73CS_LongLat_Sub_6_Block']); longLatRun = 1;
% patID = 'P73CS';
% FFAChansOnly = 1; 

% basePath = [diskPath filesep 'Face_Screening' filesep 'P73CS' filesep '2000FaceScreen_Session_1_20210326'];
% taskStruct = load([basePath filesep 'P73CS_FullSynthFaces_Sub_3_Block']);
% patID = 'P73CS';
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'AIC_Task' filesep 'P73CS' filesep 'Screening_Session_1_20210328'];
% taskStruct = load([basePath filesep 'P73CS_LL_AIC_Sub_6_Block']); patID = 'P73CS';
% FFAChansOnly = 1; 

% basePath = [diskPath filesep 'Object_Screening' filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328'];
% taskStruct = load([basePath filesep 'P73CS_ParamObj_Sub_4_Block']); patID = 'P73CS';
% FFAChansOnly = 0; 
% siphonOff500 = 1;

% basePath = [diskPath filesep 'AIC_Task' filesep 'P73CS' filesep 'ReScreenAIC_Session_1_20210328'];
% taskStruct = load([basePath filesep 'P73CS_AICReScreen_Sub_6_Block']);
% reScreenRun = 1; patID = 'P73CS';

% basePath = [diskPath filesep 'Face_Screening' filesep 'P73CS' filesep 'FaceViewScreening_Session_2_20210331'];
% taskStruct = load([basePath filesep 'P73CS_FWFast_Sub_6_Block']); patID = 'P73CS';
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'AIC_Task' filesep 'P74CS' filesep 'AICScreening_Session_1_20210823'];
% taskStruct = load([basePath filesep 'P74CS_AIC_Sub_6_Block']); patID = 'P74CS';

% basePath = [diskPath filesep 'Object_Screening' filesep 'P75CS' filesep 'FingerprintScreening_Session_1_20210902'];
% taskStruct = load([basePath filesep 'P75CS_ObjScreen_Sub_4_Block']); patID = 'P75CS';
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P76CS' filesep 'FingerprintScreening_Session_1_20210915'];
% taskStruct = load([basePath filesep 'P76CSFast_Sub_4_Block']); patID = 'P76CS';
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_1_20210917'];
% taskStruct = load([basePath filesep 'P76CSFast_2_Sub_4_Block']); patID = 'P76CS';
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];
% taskStruct = load([basePath filesep 'P76CSRec_ReScreen_Sub_4_Block']); patID = 'P76CS';
% log = [basePath filesep 'P76CSRec_ReScreen_Screening_2021-09-17_16-42-19.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925'];
% taskStruct = load([basePath filesep 'P76CS_RecScreen_2_Sub_4_Block']); patID = 'P76CS';
% log = [basePath filesep 'P76CS_RecScreen_2_Screening_2021-09-25_13-20-38.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_3_20210927'];
% taskStruct = load([basePath filesep 'P76CS_RecScreen_3_Sub_4_Block']); patID = 'P76CS';
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'];
% taskStruct = load([basePath filesep 'P76CSRec_ReScreen_3_Sub_4_Block']); patID = 'P76CS';
% % taskStruct.subID = 'P76CSRec_ReScreen_3'; % TEMPORARY
% log = [basePath filesep 'P76CSRec_ReScreen_3_Screening_2021-09-27_17-12-16.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P77CS' filesep 'FingerprintScreening_Session_1_20211007'];
% taskStruct = load([basePath filesep 'P77CS_1_Sub_4_Block']); patID = 'P77CS';
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P78CS' filesep 'FingerprintScreening_Session_1_20220309'];
% taskStruct = load([basePath filesep 'P78CS_1_Sub_4_Block']); patID = 'P78CS';
% log = [basePath filesep 'P78_Screen1_Screening_2022-03-09_15-00-56.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P78CS' filesep 'FingerprintScreening_Session_2_20220314'];
% taskStruct = load([basePath filesep 'P78CS_Screen2_Sub_4_Block']); patID = 'P78CS';
% log = [basePath filesep 'P78CS_Screen2_Screening_2022-03-14_10-36-08.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'RecallScreening_Session_1_20220330'];
% taskStruct = load([basePath filesep 'P79CS_1_Sub_4_Block']); patID = 'P79CS';
% log = [basePath filesep 'P79CS_1_Screening_2022-03-30_10-40-25.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'];
% taskStruct = load([basePath filesep 'P79CS_ReScreen_1_Sub_4_Block']); patID = 'P79CS';
% log = [basePath filesep 'P79CS_ReScreen_1_Screening_2022-03-30_15-46-44.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P79CS' filesep 'FingerprintScreening_Session_1_20220331'];
% taskStruct = load([basePath filesep 'P79CS_2_Sub_4_Block']); patID = 'P79CS';
% log = [basePath filesep 'P79CS_2_Screening_2022-03-31_11-35-11.txt'];
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'RecallScreening_Session_2_20220403'];
% taskStruct = load([basePath filesep 'P79CS_3_Sub_4_Block']); patID = 'P79CS';
% log = [basePath filesep 'P79CS_3_Screening_2022-04-03_09-24-54.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'];
% taskStruct = load([basePath filesep 'P79CS_ReScreen_3_Sub_4_Block']); patID = 'P79CS';
% log = [basePath filesep 'P79CS_ReScreen_3_Screening_2022-04-03_13-03-07.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'RecallScreening_Session_3_20220405'];
% taskStruct = load([basePath filesep 'P79CS_4_Sub_4_Block']); patID = 'P79CS';
% log = [basePath filesep 'P79CS_4_Screening_2022-04-05_17-02-37.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'];
% taskStruct = load([basePath filesep 'P79CS_ReScreen_4_Sub_4_Block']); patID = 'P79CS';
% log = [basePath filesep 'P79CS_ReScreen_4_Screening_2022-04-05_19-39-04.txt'];
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P80CS' filesep 'FingerprintScreening_Session_1_20220727'];
% taskStruct = load([basePath filesep 'P80CS_2_Sub_4_Block']); patID = 'P80CS';
% log = [basePath filesep 'P80CS_2_Screening_2022-07-27_15-35-04.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'RecallScreening_Session_1_20220728'];
% taskStruct = load([basePath filesep 'P80CS_RecScreen_1_Sub_4_Block']); patID = 'P80CS';
% log = [basePath filesep 'P80CS_RecScreen_1_Screening_2022-07-28_12-56-01.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'];
% taskStruct = load([basePath filesep 'P80CS_ReScreenRecall_Sub_4_Block']); patID = 'P80CS';
% log = [basePath filesep 'P80CS_ReScreenRecall_Screening_2022-07-28_16-53-41.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'RecallScreening_Session_2_20220730'];
% taskStruct = load([basePath filesep 'P80CS_RecScreen_2_Sub_4_Block']); patID = 'P80CS';
% log = [basePath filesep 'P80CS_RecScreen_2_Screening_2022-07-30_19-45-00.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731'];
% taskStruct = load([basePath filesep 'P80CS_ReScreecRecall_2_Sub_4_Block']); patID = 'P80CS';
% log = [basePath filesep 'P80CS_ReScreecRecall_2_Screening_2022-07-31_15-21-19.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P80CS' filesep 'FingerprintScreening_Session_2_20220801'];
% taskStruct = load([basePath filesep 'P80CS_2_Att2_Sub_4_Block']); patID = 'P80CS';
% log = [basePath filesep 'P80CS_2_Att2_Screening_2022-08-01_21-01-52.txt'];
% FFAChansOnly = 0;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'FingerprintScreening_Session_1_20221026'];
% taskStruct = load([basePath filesep '81CS_forReal_Sub_4_Block']); patID = 'P81CS';
% log = [basePath filesep '81CS_forReal_Screening_2022-10-26_10-37-01.txt'];
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'FingerprintScreening_Session_2_20221028'];
% taskStruct = load([basePath filesep 'P81CS_2_Sub_4_Block']); patID = 'P81CS';
% log = [basePath filesep 'P81CS_2_Screening_2022-10-28_16-39-34.txt'];
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030'];
% taskStruct = load([basePath filesep 'P81CS_AM_Sub_4_Block']); patID = 'P81CS';
% log = [basePath filesep 'P81CS_AM_Screening_2022-10-30_09-32-21.txt'];
% FFAChansOnly = 1;

basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030'];
taskStruct = load([basePath filesep 'P81_synth_Sub_4_Block']); patID = 'P81CS';
log = [basePath filesep 'P81_synth_Screening_2022-10-30_16-08-03.txt'];
FFAChansOnly = 1;
CLSynthScreen = 1; noNatSort = 1;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P82CS' filesep 'FingerprintScreening_Session_1_20230111'];
% taskStruct = load([basePath filesep 'P82CS_1_Sub_4_Block']); patID = 'P82CS';
% log = [basePath filesep 'P82CS_1_Screening_2023-01-11_11-22-36.txt'];
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115'];
% taskStruct = load([basePath filesep 'P82CS_CL_1_Sub_4_Block']); patID = 'P82CS';
% log = [basePath filesep 'P82CS_CL_1_Screening_2023-01-15_13-01-44.txt'];
% FFAChansOnly = 1;

% basePath = [diskPath filesep 'Object_Screening' filesep 'P82CS' filesep 'ClosedLoopReScreen_Session_1_20230115'];
% taskStruct = load([basePath filesep 'P82CS_CLReScreen_Sub_4_Block']); patID = 'P82CS';
% log = [basePath filesep 'P82CS_CLReScreen_Screening_2023-01-15_18-25-22.txt'];
% FFAChansOnly = 1;
% CLSynthScreen = 1; noNatSort = 0;



sortPath = [filesep 'sort'];
finalPath = [filesep 'final'];
rawPath = [filesep 'raw'];

% if strcmp(host(1:end-1), 'DWA644201')
%     pathTaskCode = 'D:\Users\wadiav\Dropbox\Caltech\Thesis\Human_work\Cedars\epicScreeningVarun';
% elseif strcmp(host(1:end-1), 'DESKTOP-LJHLIED')
%     pathTaskCode = 'E:\Dropbox\Caltech\Thesis\Human_work\Cedars\epicScreeningVarun';
% else % mac
%     pathTaskCode = '/Users/varunwadia/Dropbox/Caltech/Thesis/Human_work/Cedars/epicScreeningVarun';
% end
pathTaskCode = [boxPath filesep 'epicScreeningVarun'];

% basePath = [diskPath filesep 'Stories_Task' filesep 'forJie'];
% taskStruct = load([basePath filesep 'P65CS_Sub_SStories_02-02-2020_16-07-11.mat']); patID = 'P65CS';
% pathTaskCode = '/Users/varunwadia/Dropbox/Caltech/Thesis/Human_work/Cedars/storiesTaskVarun';
addpath(pathTaskCode);
setTTLCodes;

% For Deep Learning stuff
% --------------------------------------------------------------------
% pathStimuli = 'Z:\LabUsers\vwadia\screeningTaskVarun\AICScreeningStimuli_png';
% pathStimuli = 'Z:\LabUsers\vwadia\screeningTaskVarun\ResizedAICScreeningStimuli';
% pathStimuli = 'Z:\LabUsers\vwadia\screeningTaskVarun\ResizedAICScreeningStimuli';
pathStimuli = [basePath filesep 'stimuliUsed'];
% --------------------------------------------------------------------

% collecting events
% -------------------------------------------------------------------------
events = getRawTTLs([basePath filesep rawPath filesep 'Events.nev'], 1);
% if there is more than 1 start TTL
% happens if session needed to be aborted and restarted in the pt room
% if (length(find(events(:, 2) == EXPERIMENT_ON)) > 1)
%     disp('Make sure correct event indices are used!');
%     keyboard;
% end
% events = events(146:end, :); % for P69 AICScreen session 1
% events = events(73:1495, :); % for P70 AICScreen session 1
% events = events(5:988, :); % for P70 AICScreen session 2
% events = events(14:1026, :); % for P70 AICScreen session 2
% -------------------------------------------------------------------------

% if ~strcmp(taskStruct.subID, '62')
%     if strcmp(taskStruct.subID, 'P77CS_1')
%         FFAChans = [209:224];
%     elseif strcmp(taskStruct.subID, 'P78_Screen1') || strcmp(taskStruct.subID, 'P78CS_Screen2')
%         FFAChans = [25:32];
%     elseif strcmp(taskStruct.subID, 'P76CSFast_2')
%         FFAChans = [209:216];
%     end
%     
% else
%     FFAChans = [169:176]; 
% end

[strctCells, dupCells] = Utilities.extractCells([basePath filesep sortPath filesep finalPath] , basePath, [], patID, taskStruct.subID);

% forJie = [13 22 27 60 81];
% ind = logical(zeros(length(strctCells), 1));
% ind(forJie) = true;
% strctCells = strctCells(ind);
%% carve up regions if desired

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

% IT by itself
% IT+MTL
% IT+MTL+MFC
% IT+MFC
if FFAChansOnly
    IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
    strctCells = strctCells(IT_Cells);
end
% IT_MTL_Cells = cellfun(@(x) ismember(x, {'LA', 'LH', 'RA', 'RH', 'RFFA'}), strctCELL(:, 4), 'UniformOutput', false);
% strctCells = strctCells(cell2mat(IT_MTL_Cells));

% IT_MTL_MFC_Cells = cellfun(@(x) strcmp(x, {'LA', 'LH', 'RA', 'RH', 'LSMA', 'RSMA', 'RFFA'}), strctCELL(:, 4));
% strctCells = strctCells(cell2mat(IT_MTL_MFC_Cells));

% IT_MFC_Cells = cellfun(@(x) strcmp(x, ['LSMA', 'RSMA', 'RFFA']), strctCELL(:, 4));
% strctCells = strctCells(cell2mat(IT_MFC_Cells));


%% make screeningData struct

screeningData = struct;
screeningData.subID = taskStruct.subID;
screeningData.psth = cell(length(strctCells), 3);
screeningData.responses = cell(length(strctCells), 1);
% screeningData.magnitudeOrder = cell(length(strctCells), 1); % order of images sorted by response magnitude

% find(events(:, 2) == 0) % had to disconnect the pci-e cable again
% spits out 2, 6, 768, 769 for P63 Screening 2
experimentBegin = find(events(:, 2) == EXPERIMENT_ON);
experimentEnd = find(events(:, 2) == EXPERIMENT_OFF);

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
    %         names = cell2mat(names);
    %
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

% if strcmp(taskStruct.subID, '81CS_forReal')
%     IMAGE_ON = IMAGE_ON_LOG; IMAGE_OFF = IMAGE_OFF_LOG;
% end

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

%% extract rasters and psths

% for raster extraction and plotting
% ------------------------------------------------------------------------
if stimOffDur < stimDur
    tl_1 = round((stimOffDur)*1e-3, 2) %'tee el dash one'
else
    tl_1 = round((stimDur)*1e-3, 2) %'tee el dash one'
end
if tl_1 == 0.13
    tl_1 = 0.17
end
tl_2 = round((stimDur)*2e-3, 2)
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
  
%     t_rem = [];
%     for i = 1:length(ord500Img)
%         t = find(sorted500Img == i);
%         if length(t) > screeningData.numRepetitions
%             t_rem = [t_rem; t(end)]; % should I remove the last rep or some random one in the middle?
%         end
%     end
%     correct500Img(t_rem) = []; % removing repetitions
%     
%     % so that anova labels are made correctly
%     screeningData.imageIDs = [1:500]';
%     screeningData.sortedOrder = ord500Img(correct500Img); 
% %     assert(isequal(screeningData.sortedOrder, repelem(screeningData.imageIDs, screeningData.numRepetitions)), 'IMAGE ORDER IS WRONG');
    
end

%% siphoning off synthetic images from the closed loop screen

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


%% ANOVA for selectivity - make the labels in a function

% making for different anova types
% -----------------------------standard------------------------------------
% anovaType = '1xN';
% labels = sortedOrder;
% -------------------------------Category----------------------------------
[labels, anovaType] = Utilities.makeObjCatLabsScreening(taskStruct.subID, screeningData.sortedOrder);
labels(labels == 0) = [];

%% create labels for quadrant specificity anova

% %       |
% %    2  |  1
% % ------+------
% %    3  |  4
% %       |
% load([diskPath filesep 'ObjectSpace' filesep '1593Stimuli' filesep 'params_Alexnet_fc6_1593_objects.mat']);
% load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']);
% if exist('params')
%     score = params;
% end
% PC12 = score(:, 1:2);
% Q1 = []; Q2 = []; Q3 = []; Q4 = [];
% for ii = 1:length(PC12(:, 1))
%     if (PC12(ii, 1) > 0 && PC12(ii, 2) > 0)
%         Q1 = vertcat(Q1, ii);
%     elseif (PC12(ii, 1) < 0) && (PC12(ii, 2) > 0)
%         Q2 = vertcat(Q2, ii);
%     elseif (PC12(ii, 1) < 0) && (PC12(ii, 2) < 0)
%         Q3 = vertcat(Q3, ii);
%     else (PC12(ii, 1) > 0) && (PC12(ii, 2) < 0);
%         Q4 = vertcat(Q4, ii);
%     end
% end
% 
% 
% % labels
% PC12(Q1, 3) = 1;
% PC12(Q2, 3) = 2;
% PC12(Q3, 3) = 3;
% PC12(Q4, 3) = 4;
% 
% % this label assignment works even with 1-back extras
% catOrder = zeros(length(screeningData.sortedOrder), 1);
% for i = 1:length(screeningData.imageIDs)
%     catOrder(find(screeningData.sortedOrder == screeningData.imageIDs(i))) = PC12(i, 3);
% end
% labels = catOrder;
% anovaType = 'CategoryQuadrant';

%% compute response latency for each cell and responses to each image 

% screeningData.responses = cell(length(strctCells), 3);

% computing response latency and calculating responses
if exist('catOrder')
    labels = catOrder;
end

% which to use 
% 0 - basic, 
% 1 - sliding window anova, 
% 2 - peak of omega squared (also sliding window anova)
% 3 - Poisson spike train analysis
method = 1; 
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
                screeningData.sortedOrder, screeningData.timelimits, screeningData.stimDur, true);
            resplat = -screeningData.timelimits(1)*1e3 + respLat; % adjust to be of the same for mas the other methods
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
        
        % for Hristos
        %         trial_resp{cellIndex, 1} = sum(screeningData.psth{cellIndex, 1}(:, windowBegin:windowEnd), 2);
        
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
% responses = screeningData.responses;
% psths = screeningData.psth;
% order = screeningData.sortedOrder;

% index = cellfun(@isempty, responses);
% responses(index(:, 1), :) = [];
% psths(index(:, 1), :) = [];
% strctCells(index(:, 1)) = [];
% save([basePath filesep 'PsthandResponses'], 'screeningPsth', 'responses', 'order');
% save([basePath filesep 'strctCells'], 'strctCells');
% save([basePath filesep 'ResponseLatencies'], 'responses');
% save([basePath filesep 'ITCells_1593Stim_Scrn_basicMethod_2Stds'], 'responses', 'psths', 'strctCells', '-v7.3')
%% Finding stim for rivalry

% [stim1BR, stim2BR, responseMatrix, screeningData] = findStimPairBR(windowBegin, windowEnd, windowLength, screeningData);

% xcorr work
% 90, 146

%% CLScreen plotting all images

% % put synth images back in the lineup
% screeningData.psth(:, 1:2) = cellfun(@(x, y) [x;y], screeningData.psth(:, 1:2), screeningData.synthPsth(:, 1:2), 'UniformOutput', false);
% 
% anovaType = '1xN';
% labels = sortedOrder;


%% Actually compute ANOVA
if ~strcmp(anovaType, '1xN') 
    screeningData.catOrder = labels;
    screeningData.catIDs = unique(labels);
end
screeningData = Utilities.selectivityCheck_ANOVA(labels, anovaType, strctCells, screeningData);
screeningData.anovaType = anovaType;

keyboard
%%
% for paper - vwadia March 2022
responses = screeningData.responses;
screeningPsth = screeningData.psth;

if ~exist('CLSynthScreen', 'var') || CLSynthScreen == 0
    order = screeningData.sortedOrder;
else
    % only the real images
    order = screeningData.sortedOrder(1:(length(screeningData.imageIDs)*screeningData.numRepetitions));    
end

if FFAChansOnly
    save([basePath filesep 'PsthandResponses'], 'screeningPsth', 'responses', 'order', '-v7.3')
    save([basePath filesep 'strctCells'], 'strctCells')
    
elseif exist('IT_MTL_Cells')
    save([basePath filesep 'IT_MTL_PsthandResponses'], 'screeningPsth', 'responses', 'order', '-v7.3')
    save([basePath filesep 'IT_MTLCells'], 'strctCells')
    
else
    save([basePath filesep 'PsthandResponses'], 'screeningPsth', 'responses', 'order', '-v7.3')
    save([basePath filesep 'strctCells'], 'strctCells')
    
end

if exist('CLSynthScreen', 'var') && CLSynthScreen == 1
    synthResponses = screeningData.synthResponses;
    synthPsth = screeningData.synthPsth;
    synthOrder = synthOrd_postSort;
    
    save([basePath filesep 'SynthPsthandResponses'], 'synthPsth', 'synthResponses', 'synthOrder', '-v7.3')
end

% Screening.scrn_rollingwindow_STA_analysis

keyboard;

% save([basepath filesep 'dataForKevin'], 'strctCells', 'screeningData', '-v7.3')
% save([basePath filesep 'dataForJuri_Autocorr'], 'psths', 'responses', 'order',  'strctCells', '-v7.3')


%% for Hristos - STA ground truth cell

% % cell 8276 Obj screen
% psths = screeningData.psth(:, :);
% % responses = screeningData.responses{4, 1};
% psthOrder = screeningData.sortedOrder;
% 
% 
% save([basePath filesep 'P73CS_forHristos'], 'strctCells', 'psthOrder', 'trial_resp');
% 
%% for Umit
% 
% stimOrder = order;
% psthOrder = screeningData.sortedOrder;
% psth = screeningData.psth;
% 
% save([basePath filesep 'P73CS_CategoryScreen2_FullBrain_forUmit'], 'strctCells', 'psth', 'psthOrder', 'stimOrder', 'events', '-v7.3');

% categoryOrder = screeningData.catOrder;
% categoryIDs = {1, 'Faces'; 2, 'Text'; 3, 'Flora'; 4, 'Fauna'; 5, 'Objects'};
% Objfeatures = load(['G:\SUAnalysis\ObjectSpace' filesep 'parameters_1593_objects.mat']);
% features_1593_objects = Objfeatures.score;
% save([basePath filesep 'P73CS_ObjectScreen_CategoryandFeatures_forUmit'], 'categoryOrder', 'categoryIDs', 'features_1593_objects');
% 
% save([basePath filesep 'P73CS_ObjectScreen_CategoryandFeatures_forUmit'], 'parameters_2k_synthetic_faces.mat');



%% ranksum test for selectivity - check if this is legit

% offset = 50;
% windowLength = ceil(stimDur);
% windowBegin = tl_1*1e3 + offset;
% windowEnd = windowBegin+windowLength;
% [screeningData] = responsivityCheck(windowBegin, windowEnd, windowLength, offset, strctCells, screeningData, 0.05);

%% Plotting - all rasters
sortByMag = 0;


if (isfield(screeningData, 'catIDs') && ~isempty(screeningData.catIDs))
    if strcmp(anovaType, 'CategoryObject') 
        if strcmp(taskStruct.subID, 'P71CS')
            % Category
            screeningData.lgnd = {'Faces', 'Bodies', 'Fruits', 'Hands', 'Techno'};
        elseif strcmp(taskStruct.subID, 'P73CS')
            % Category
            screeningData.lgnd = {'Faces', 'Names (text)', 'Places', 'Scrambled', 'Objects'};
        elseif strcmp(taskStruct.subID, 'P73CS_FWFast')
            % Face View
            screeningData.lgnd = {'front', 'halfLeft', 'fullLeft', 'halfRight', 'fullRight', 'lookUp', 'lookDown', 'backOfHead', 'objects'};
        elseif strcmp(taskStruct.subID, 'P73CS_ParamObj') || strcmp(taskStruct.subID, 'P73CS_ParamObj_500') || strcmp(taskStruct.subID, 'P71_large') || strcmp(taskStruct.subID, 'P71CS_Fast') || strcmp(taskStruct.subID, 'P71CS_Object')...
                || strcmp(taskStruct.subID, '62') || strcmp(taskStruct.subID, 'P75CS_ObjScreen') || strcmp(taskStruct.subID, 'P76CSFast') || strcmp(taskStruct.subID, 'P76CSFast_2')...
                || strcmp(taskStruct.subID, 'P76CSRec_ReScreen') || strcmp(taskStruct.subID, 'P76CS_RecScreen3') || strcmp(taskStruct.subID, 'P76CS_RecScreen_3') || strcmp(taskStruct.subID, 'P76CSRec_ReScreen_3')...
                || strcmp(taskStruct.subID, 'P77CS_1') || strcmp(taskStruct.subID, 'P78_Screen1') || strcmp(taskStruct.subID, 'P78CS_Screen2') || strcmp(taskStruct.subID, 'P79CS_1')...
                || strcmp(taskStruct.subID, 'P79CS_2') || strcmp(taskStruct.subID, 'P79CS_ReScreen_1') || strcmp(taskStruct.subID, 'P79CS_3') || strcmp(taskStruct.subID, 'P79CS_ReScreen_3')...
                || strcmp(taskStruct.subID, 'P79CS_4') || strcmp(taskStruct.subID, 'P79CS_ReScreen_4') || strcmp(taskStruct.subID, 'P80CS_2') || strcmp(taskStruct.subID, 'P80CS_RecScreen_1')...
                || strcmp(taskStruct.subID, 'P80CS_ReScreenRecall') || strcmp(taskStruct.subID, 'P80CS_RecScreen_2') || strcmp(taskStruct.subID, 'P80CS_2_Att2') || strcmp(taskStruct.subID, '81CS_forReal')...
                || strcmp(taskStruct.subID, 'P81CS_2') || strcmp(taskStruct.subID, 'P81CS_AM') || strcmp(taskStruct.subID, 'P81_synth') || strcmp(taskStruct.subID, 'P82CS_1') || strcmp(taskStruct.subID, 'P82CS_CL_1')...
                || strcmp(taskStruct.subID, 'P82CS_CLReScreen')
            % 1593 Param Objs
            screeningData.lgnd = {'Faces', 'Text', 'Plants/fruits', 'Animals', 'Objects'};
%             screeningData.lgndCols = 
        elseif  strcmp(taskStruct.subID, 'P73CS_AICReScreen') || strcmp(taskStruct.subID, 'P74CS_AIC')
            % Re-Screen AIC
            screeningData.lgnd = {'Faces', 'Names (text)', 'Places', 'Objects'};
        elseif strcmp(taskStruct.subID, 'P71CS_RecScreen2')
            % screening with 1224 object set
            screeningData.lgnd = {'Faces', 'Plants/fruits', 'Animals', 'Objects'};            
        end
    elseif strcmp(anovaType, 'CategoryQuadrant')
        screeningData.lgnd = {'Q1-SpikyInanimate', 'Q2-StubbyInanimate', 'Q3-StubbyAnimate', 'Q4-SpikyAnimate'};
    elseif strcmp(anovaType, 'SingleCategory')
        
        %         screeningData.lgnd = {'1593 Param Objects'};%
        screeningData.lgnd = {'Synthesized Faces'};%
        
    end
    imagesPerSubplot = length(screeningData.catIDs);
    subPlotNum = 1;
    numFigsPerCell = 1;
else
    imagesPerSubplot = 10;
    subPlotNum = 10;
    numFigsPerCell = 6; %ceil(length(screeningData.imageIDs)/(imagesPerSubplot*subPlotNum)); % each figure has max 100 images plotted
end

% plot and save
Utilities.Plotting.PlotRastersScreening(basePath, imagesPerSubplot, subPlotNum, numFigsPerCell, strctCells, screeningData, sortByMag)

% keyboard

%% post it plots
onlyTopandBottom = 1;
for cellIndex = l(strctCells)
    
    if ~exist('CLSynthScreen', 'var') || CLSynthScreen == 0
        sortOrd = screeningData.sortedOrder;
    else
        % only the real images
        sortOrd = screeningData.sortedOrder(1:(length(screeningData.imageIDs)*screeningData.numRepetitions));
    end
    psth = screeningData.psth(cellIndex, :);
    timelimits = screeningData.timelimits;
    % note that respLat here has had the timelimit subtracted already
    emptyIndex = cellfun('isempty', screeningData.responses(:, 3));     % Find indices of empty cells
    screeningData.responses(emptyIndex, 3) = {0};                    % Fill empty cells with 0
    cellInQuestion = find(cell2mat(screeningData.responses(:, 3))==strctCells(cellIndex).Name);
    
    if ~isempty(cellInQuestion)
        respLat = screeningData.responses{cellInQuestion, 2};
        cellInfo = strctCells(cellIndex);
        
        if onlyTopandBottom
            options.topStimNum = 20;
            rows = 4;
            cols = 20;
            handles = Utilities.Plotting.makePostItPlot(psth, sortOrd, timelimits, respLat-screeningData.timelimits(1)*1e3, rows, cols, basePath, pathStimuli, cellInfo, screeningData, options);
        else
            rows = 10;
            cols = 20;
            handles = Utilities.Plotting.makePostItPlot(psth, sortOrd, timelimits, respLat-screeningData.timelimits(1)*1e3, rows, cols, basePath, pathStimuli, cellInfo, screeningData);
        end
        
    end
    
end

%% STA plots

if strcmp(taskStruct.subID, 'P73CS_Full') 
%     z_scored = 0;
    load([diskPath filesep 'ObjectSpace' filesep 'parameters_2k_synthetic_faces.mat']);
    params = params(1:667, :);
    options.screenType = 'Face';
elseif strcmp(taskStruct.subID, 'P73CS_ParamObj') || strcmp(taskStruct.subID, '62')
    load([diskPath filesep 'ObjectSpace' filesep '1593Stimuli' filesep 'params_Alexnet_fc6_1593Stimuli.mat']); % will create score = 1593x50
    options.screenType = 'Object';
elseif strcmp(taskStruct.subID, 'P75CS_ObjScreen') || strcmp(taskStruct.subID, 'P76CSFast') || strcmp(taskStruct.subID, 'P76CSFast_2')...
        || strcmp(taskStruct.subID, 'P76CSRec_ReScreen') || strcmp(taskStruct.subID, 'P76CS_RecScreen3') || strcmp(taskStruct.subID, 'P76CS_RecScreen_3')...
        || strcmp(taskStruct.subID, 'P76CSRec_ReScreen_3') || strcmp(taskStruct.subID, 'P77CS_1') || strcmp(taskStruct.subID, 'P78_Screen1')...
        ||  strcmp(taskStruct.subID, 'P78CS_Screen2') || strcmp(taskStruct.subID, 'P79CS_1') || strcmp(taskStruct.subID, 'P79CS_2') || strcmp(taskStruct.subID, 'P79CS_ReScreen_1')...
        || strcmp(taskStruct.subID, 'P79CS_3') || strcmp(taskStruct.subID, 'P79CS_ReScreen_3') || strcmp(taskStruct.subID, 'P79CS_4') || strcmp(taskStruct.subID, 'P79CS_ReScreen_4')...
        || strcmp(taskStruct.subID, 'P80CS_2') || strcmp(taskStruct.subID, 'P80CS_RecScreen_1') || strcmp(taskStruct.subID, 'P80CS_ReScreenRecall') || strcmp(taskStruct.subID, 'P80CS_RecScreen_2')...
        || strcmp(taskStruct.subID, 'P80CS_2_Att2') || strcmp(taskStruct.subID, '81CS_forReal') || strcmp(taskStruct.subID, 'P81CS_2') || strcmp(taskStruct.subID, 'P81CS_AM') || strcmp(taskStruct.subID, 'P81_synth')...
        || strcmp(taskStruct.subID, 'P82CS_1') || strcmp(taskStruct.subID, 'P82CS_CL_1')
    
        load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create params = 500x50
        options.screenType = 'Object';
elseif strcmp(taskStruct.subID, 'P73CS_ParamObj_500')
    % get params matrix
%     params2 = Utilities.getDNActivations(path500Stimuli, screeningData.imageIDs, 50, 'fc6');
    options.screenType = 'Object';
    load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create params = 500x50
%     assert(isequal(params, params2), 'WRONG IMAGES');
else
    % get params matrix
    params = Utilities.getDNActivations(pathStimuli, screeningData.imageIDs, 50, 'fc6');
    options.screenType = 'Object';
    disp('hi');
end
% screeningData.imageIDs = 1:500;
% basePath = taskPath;
% screeningData.responses = responses;
% pathOut = [basePath filesep 'STA_and_projections_stdNorm'];
pathOut = [basePath filesep 'STA_and_projections_varWindow'];

if ~exist(pathOut)
    mkdir([pathOut]);
end

% % Cell 4 - P73CS param Obj screening 1593 stimuli
% options.encoded_stim = [242 298 900];
% options.encoded_stim_ortho = [403 793]
% % Cell 4 - P73CS param Obj screening 500 stimuli
% options.encoded_stim = [463 425 403]; %[6 235 497] xy positions
% options.encoded_stim_ortho = [427 238]; %[286 64] positions
% ---------------------------------------------
% % Cell 22 - P73CS param Obj screening 1593 stimuli
% options.encoded_stim = [143 439 546];
% options.encoded_stim_ortho = [246 1069];
% % Cell 22 - P73CS param Obj screening 500 stimuli
% options.encoded_stim = [408 207 170]; %[157 434 497] xy positions
% options.encoded_stim_ortho = [434 263]; %[162 117]
% -----------------------------------------
% % Cell 33 - P73CS param face screen
% options.encoded_stim = [639 298 45];
% options.encoded_stim_ortho = [223 600];

% options.encoded_stim = [210 405 414];
% options.encoded_stim = [489 431 226 330 16 218];
if strcmp(taskStruct.subID, 'P76CSRec_ReScreen') || strcmp(taskStruct.subID, 'P76CSFast_2')
    options.marked_positions = [1 456 457 221 376 114]; % P76 Recall 1
    options.recalled_stim = [12 19 25 123 270 487]; % P76 Recall 1
elseif strcmp(taskStruct.subID, 'P76CS_RecScreen3')
    options.marked_positions = [16 34 144 358 382 450]; % P76 Recall 2
    options.recalled_stim = [54 129 130 186 270 449]; % P76 Recall 2
elseif  strcmp(taskStruct.subID, 'P76CSRec_ReScreen_3') || strcmp(taskStruct.subID, 'P76CS_RecScreen_3')
    options.marked_positions = [175 107 345 340 476 459 499 496]; % P76 Recall 3
    % options.recalled_stim = [230 344 81 45 135 181 44 18]; % P76 Recall 3
    options.recalled_stim = [18 44 45 81 135 181 230 344];
elseif  strcmp(taskStruct.subID, 'P79CS_1') || strcmp(taskStruct.subID, 'P79CS_ReScreen_1')
    options.marked_positions = [3 104 167 196 422 453 473 491];
    options.recalled_stim = [9 157 167 200 201 291 422 498];
elseif  strcmp(taskStruct.subID, 'P79CS_3') || strcmp(taskStruct.subID, 'P79CS_ReScreen_3')
    options.marked_positions = [476 497 453 495 279 309 13 167];
    options.recalled_stim = [9 12 117 292 360 368 421 492];
elseif  strcmp(taskStruct.subID, 'P79CS_4') || strcmp(taskStruct.subID, 'P79CS_ReScreen_4')
    options.marked_positions = [49 60 70 173 422 480 499 500];
    options.recalled_stim = [77 112 160 232 278 345 387 440];
elseif  strcmp(taskStruct.subID, 'P80CS_RecScreen_1') || strcmp(taskStruct.subID, 'P80CS_ReScreenRecall')
    options.marked_positions = [ 5   185   388   431   462   468   486   498];
    options.recalled_stim = [17    61    76   114   157   161   177   480 ];
elseif strcmp(taskStruct.subID, 'P80CS_RecScreen_2') || strcmp(taskStruct.subID, 'P80CS_ReScreecRecall_2') 
    options.marked_positions = [134 257 454 466 488 492 498 499];
    options.recalled_stim = [55 88 148 251 256 274 285 365];
end


% write a script to check images automatically (load in stimUsedRecall and
% general imDir and compare filenames after sampling recalled_stim from
% imDir

% 77   112   160   232   278   345   387   440
% options.encoded_stim = [134:210];


options.recalledCols = [1 0.5 0;...% orange
    1 0.25 0;...% orange
    0.8 1 0;...% yellow
    0 1 0;...% green
    0 1 1;...% turquoise
    1 0 1;...% pink
    0.75 0.75 0.75;...% grey
    0 0 0];% white 
% options.recalledCols = [1 0.5 0;...% orange
%     0.8 1 0;...% yellow
%     0 1 0;...% green
%     0 1 1;...% turquoise
%     1 0 1;...% pink
%     0.75 0.75 0.75];% grey 

alpha = 0.01;
for cellIndex = l(strctCells)
    
    if ~isempty(screeningData.responses{cellIndex, 2})
        options.ind_train = screeningData.imageIDs; % use all objects to calculate STA
%         options.encoded_stim = [19, 21, 69, 3, 51]; % caharacters in AIC screen March28th
%         [hfig, maxCoords, minCoords, medCoords, max_stim, min_stim, med_stim] = STA_figure(screeningData.responses{cellIndex}, score, options, z_scored); % pass score to this instead of projectedResponses
        [hfig, p, options] = Utilities.ObjectSpace.STA_figure_original(screeningData.responses{cellIndex}, params, options); % pass score to this instead of projectedResponses
        if isfield(options, 'xvals') && isfield(options, 'yvals')
            strctCells(cellIndex).Im_xvals = options.xvals;
            strctCells(cellIndex).Im_yvals = options.yvals;
        end
        if isfield(options, 'recalled_stim')
            strctCells(cellIndex).recalledStim = options.recalled_stim;
        end
%         [hfig] = face_id_analysis_STA_example1(screeningData.responses{cellIndex}, score, options); % pass score to this instead of projectedResponses
        sgtitle({['Cell number ' num2str(strctCells(cellIndex).Name)] 'STA and projections ', strctCells(cellIndex).brainArea});
        
        if p < alpha
            newPathOut = [pathOut filesep 'significant_cells'];
            if ~exist(newPathOut)
                mkdir([newPathOut]);
            end
            if isfield(options, 'encoded_stim') || isfield(options, 'recalled_stim')
                print(hfig, [newPathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_orderedStim'], '-dpng', '-r0')
            else
                print(hfig, [newPathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name)], '-dpng', '-r0')
            end
        else
            print(hfig, [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name)], '-dpng', '-r0')
        end

        %     strctCells(cellIndex).maxCoords = maxCoords;
        %     strctCells(cellIndex).minCoords = minCoords;
        %     strctCells(cellIndex).medCoords = medCoords;
        %     strctCells(cellIndex).maxStimImages = screeningData.imageIDs(max_stim);
        %     strctCells(cellIndex).minStimImages = screeningData.imageIDs(min_stim);
        %     strctCells(cellIndex).medStimImages = screeningData.imageIDs(med_stim);
        %     keyboard
        close all;
    end
end


%% Recipe for marking recall stimuli
%% Sept 25th, 2021

% updated Nov 16th, 2021
% Make plot and click on points to use (breakpoint in STA_figure_original - right after scatter is made ~line 157)
% Use those XY coordinates to find out what indicies in x, y they correspond to - see below for converting cursor info into pot_rec_stim
% in STA_figure_original you now have pot_stim_pos = options.marked_positions
% reorder_ind(indices) are the actual stimuli numbers

%% March 22nd, 2022

% As stated above, make plot and click on points to use (breakpoint in STA_figure_original - right after scatter is made ~line 157)
% Save each pair as 'cursor_info_pairx' 
% Once this is done for all 3/4 pairs dbquit and return here - use the cell below to save them as pot_rec_stim
% SAVE THE VARIABLE options.pot_rec_stim
% Re-run the STA cell above - with a breakpoint in STA_figure_original after the if statement 'exist(pot_rec_stim)'
% Note down the 'pot_stim_pos' you will get - reorder_ind(pot_stim_pos) are the actual stimuli numbers - add them to the STA cell above

%%

% for i = 1:length(cursor_info)
%     options.pot_rec_stim(i, :) = cursor_info(i).Position;
% end
for i = 1:length(cursor_info_pair1)
    options.pot_rec_stim(i, :) = cursor_info_pair1(i).Position;
end
for i = 1:length(cursor_info_pair2)
    
    options.pot_rec_stim(i+2, :) = cursor_info_pair2(i).Position;
end
for i = 1:length(cursor_info_pair3)
    options.pot_rec_stim(i+4, :) = cursor_info_pair3(i).Position;
end
for i = 1:length(cursor_info_pair4)
    options.pot_rec_stim(i+6, :) = cursor_info_pair4(i).Position;
end


% important step
pot_rec_stim = options.pot_rec_stim;
save([basePath filesep 'markedPositionsRecallStim'], 'pot_rec_stim');

%%

for i = 1:length(cursor_info_testpair)
    options.test(i, :) = cursor_info_testpair(i).Position;
end

%% stim for IT paper

options.encoded_pos(1, :) = cursor_info_o1.Position
options.encoded_pos(2, :) = cursor_info_o2.Position
options.encoded_pos(3, :) = cursor_info_p1.Position
options.encoded_pos(4, :) = cursor_info_p2.Position
options.encoded_pos(5, :) = cursor_info_p3.Position

options.pot_rec_stim = options.encoded_pos;

%%

options.encoded_pos(1, :) = cursor_info_t1.Position
options.encoded_pos(2, :) = cursor_info_t2.Position
options.pot_rec_stim = options.encoded_pos;


%% checking for number of dimensions the cell is significantly tuned to
% fr1 = screeningData.responses{4}; % example cell 8276 object screening
% for i = 1:size(score,2)
% [cc pp] = corrcoef(score(:,i),fr1);
% c(i) = cc(1,2);p(i) = pp(1,2);
% end


%%
%{
%% finding the best stim pair for BR (option 2) - images that modulated most and least cells

[so, co] = sortrows(sigImagesPop');
sigImageCount = [];
for i = unique(sigImagesPop)
    % now sigimageCount(i) contains the number of cells in which image(i)
    % evoked a sig response
    sigImageCount(end+1) = length(find(so == i));
end
stim1BR = find(sigImageCount == max(sigImageCount));
stim2BR = find(sigImageCount == min(sigImageCount));
% assert(length(stim1BR) == 1); assert(length(stim2BR) == 1);


for j = l(strctCells)
    images = screeningData.sigImagesPerCell{j};
    for k = images
        screeningData.sigCellsPerImage{k, 1}(end+1) = j;
        screeningData.sigCellsPerImage{k, 2}{end+1} = strctCells(j).brainArea;
    end
end

%% finding only significant images


screeningData.sigImagesPerCell = cell(length(strctCells), 2);
screeningData.sigCellsPerImage = cell(length(screeningData.imageIDs), 2);
sigImagesPop = [];
repetitions = 6;


% doing stats
% t-tests to see if it is visually responsive
% per cell
for cellIndex = l(strctCells)
    % per stimulus
    for imageIndex = 1:length(screeningData.imageIDs)
        
        stimOFF = screeningData.psth{cellIndex, 2}(((abs(screeningData.timelimits(1))*1000)-round(stimOffDur)):(abs(screeningData.timelimits(1))*1000),...
            ((repetitions*(imageIndex-1))+1):((repetitions*(imageIndex-1))+repetitions));
        
        stimON = screeningData.psth{cellIndex, 2}(501:1200,...
            ((repetitions*(imageIndex-1))+1):((repetitions*(imageIndex-1))+repetitions));
        [rejectNull, p] = ttest(mean(stimOFF, 1), mean(stimON, 1)); % make sure this is across trials
        
        if (rejectNull == 1 & p < 0.05)
            screeningData.sigImagesPerCell{cellIndex, 1}(end+1) = imageIndex; % we want to plot this image
            screeningData.sigImagesPerCell{cellIndex, 2}(end+1) = p;
            sigImagesPop(end+1) = imageIndex;
            
        end
    end
    
end
%}

%% plotting specific images only - colors chosen to match desired stim

% setting up viewing parameters
MarkerSize = 4;
Fontsize = 9;
numFigsPerCell = 1;
subPlotNum = 1;
imagesPerSubplot = 2;

cellIndex = 11;
imgs = [236 405];
colors = Utilities.distinguishable_colors(length(imgs));

specificOrder = [];
for i = 1:length(imgs)
    specificOrder = [specificOrder; find(screeningData.sortedOrder == imgs(i))]; 
end
specificRaster = screeningData.psth{cellIndex, 1}(specificOrder, :);
specificPSTH = screeningData.psth{cellIndex, 2}(specificOrder, :);

newSpecificOrder = repelem(1:length(imgs), screeningData.numRepetitions);
% find global ylim
% totalTrials = length(imgs)*screeningData.numRepetitions;%numFigsPerCell*imagesPerSubplot*subPlotNum;
% globalyl = Utilities.Plotting.findingGlobalYLim(specificPSTH, imgs, repelem(1:length(imgs), screeningData.numRepetitions), 'Screening', totalTrials);

for figNum = 1:numFigsPerCell
    f = figure;
    %     set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
    clf reset
    
    if iscell(strctCells(cellIndex).brainArea)
        sgtitle({[num2str(strctCells(cellIndex).Name) ' ' char(strctCells(cellIndex).brainArea) '\_' num2str(figNum)]}); % backslash allows you to print the underscore
    else
        sgtitle({[num2str(strctCells(cellIndex).Name) ' ' strctCells(cellIndex).brainArea '\_' num2str(figNum)]}); % backslash allows you to print the underscore
    end
    pathOut = [basePath filesep 'rasters'];
    if ~exist(pathOut)
        mkdir([pathOut]);
    end
    
    
    
    % currently this is redrawing the same few images on each figure, change it
    for ctr = 1:subPlotNum % make sure this is divisible - chunks of images
        
        % adjusted for which figure it is eg. if I'm ploting 100 images
        % per figure then fig 2 should start plotting from 101-200
        figurePlotOffset = ((figNum-1)*imagesPerSubplot*subPlotNum);
        imagesTOplot = 1:length(imgs);
        %         imagesTOplot=screeningData.imageIDs(((imagesPerSubplot*(ctr-1))+figurePlotOffset+1):...
        %             screeningData.imageIDs((imagesPerSubplot*(ctr-1))+figurePlotOffset+imagesPerSubplot));
        screeningData.imageIDstoDEL=setdiff(screeningData.imageIDs,[imagesTOplot]);
        
        % FR
        h_1(ctr) = subplot(2, subPlotNum, ctr);
        hold on
        
        
        for p1 = l(imagesTOplot)
            Utilities.stdshade5(specificPSTH(find(newSpecificOrder == imagesTOplot(p1)), :), 0.1, colors(mod(p1, length(imagesTOplot))+1, :), times, 2);
        end
        %         xlim([screeningData.timelimits(1)-(screeningData.timelimits(1)*0.1) screeningData.timelimits(2)-(screeningData.timelimits(2)*0.1)]);
        xlim([-0.2 0.8])
        % needs to be improved
%         ylim([0 globalyl+2]);
        yl = ylim;
        plot([0 0], [0 yl(2)+2], '--k', 'LineWidth', 1);
        plot([screeningData.stimDur/1000 screeningData.stimDur/1000], [0 yl(2)+2], '--k', 'LineWidth', 1);
        
        % raster
        h_2(ctr) = subplot(2, subPlotNum, ctr+subPlotNum);
        hold on
        for p2 = l(imagesTOplot)
            iter = find(screeningData.sortedOrder == imagesTOplot(p2)); % iter is now a vector of indices size 6
            for k = 1:screeningData.numRepetitions
                try
                    % if sorted raster used replace original with presentations
                    plot((find(specificRaster(iter(k), :)==1).*(1/1000)+screeningData.timelimits(1)),...
                        screeningData.numRepetitions*(imagesTOplot(p2)-1)+k,'Marker','square', 'LineStyle','none','MarkerFaceColor',colors(mod(imagesTOplot(p2), length(imagesTOplot))+1, :),...
                        'MarkerEdgeColor','none','MarkerSize',MarkerSize)
                    
                    hold on
                    
                end
            end
        end
        
        
        set(gca,'XGrid','on')
        set(gca, 'YGrid','on')
        grid on
        
        tick=round(length(screeningData.correctOrder)*0.20);
        
        set(gca,'FontSize',Fontsize)
        
        %         xlim([screeningData.timelimits(1)-(screeningData.timelimits(1)*0.1) screeningData.timelimits(2)-(screeningData.timelimits(2)*0.1)]);
        xlim([-0.2 0.8])
        
        ylabel('trial nr (re-ordered)','FontSize',Fontsize);
        ylim([screeningData.numRepetitions*(imagesTOplot(1)-1) screeningData.numRepetitions*(imagesTOplot(end))+1]); % why does this screw up the ratio???
        yl = ylim;
        plot([0 0], [yl(1) yl(2)], '--k', 'LineWidth', 1);
        plot([screeningData.stimDur/1000 screeningData.stimDur/1000], [yl(1) yl(2)], '--k', 'LineWidth', 1);
        xlabel('time (sec)','FontSize',Fontsize);
        
        set(gca,'FontSize',Fontsize)
        %             keyboard;
        
        
    end
    linkaxes(h_1, 'y');
    
    % saving
    filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_specificStim'];
    if ~strcmp(class(filename), 'cell')
        print(f,filename ,'-dpng','-r0')
    else
        filename = [pathOut filesep strctCells(cellIndex).brainArea{1} '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_' num2str(figNum)];
        print(f,filename ,'-dpng','-r0')
    end
    %     keyboard;
end
%}
%% Checking the PC1 and PC2 spread of preferred and least preferred stimuli

if strcmp(taskStruct.subID, 'P73CS_Full') 
%     z_scored = 0;
    load([diskPath filesep 'ObjectSpace' filesep 'parameters_2k_synthetic_faces.mat']);
    params = params(1:667, :);
    options.screenType = 'Face';
elseif strcmp(taskStruct.subID, 'P73CS_ParamObj') || strcmp(taskStruct.subID, '62')
    load([diskPath filesep 'ObjectSpace' filesep 'parameters_1593_objects.mat']); % will create score = 1593x50
    params = score;
    options.screenType = 'Object';
elseif strcmp(taskStruct.subID, 'P75CS_ObjScreen') || strcmp(taskStruct.subID, 'P76CSFast') || strcmp(taskStruct.subID, 'P76CSFast_2')...
        || strcmp(taskStruct.subID, 'P76CSRec_ReScreen') || strcmp(taskStruct.subID, 'P76CS_RecScreen3') || strcmp(taskStruct.subID, 'P76CS_RecScreen_3')...
        || strcmp(taskStruct.subID, 'P76CSRec_ReScreen_3') || strcmp(taskStruct.subID, 'P77CS_1') || strcmp(taskStruct.subID, 'P78_Screen1')...
        ||  strcmp(taskStruct.subID, 'P78CS_Screen2') ||strcmp(taskStruct.subID, 'P80CS_ReScreenRecall') || strcmp(taskStruct.subID, '81CS_forReal') 
        load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create params = 500x50
        options.screenType = 'Object';
elseif strcmp(taskStruct.subID, 'P73CS_ParamObj_500')
    % get params matrix
    params2 = Utilities.getDNActivations(path500Stimuli, screeningData.imageIDs, 50, 'fc6');
    options.screenType = 'Object';
    load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create params = 500x50
    assert(isequal(params, params2), 'WRONG IMAGES');
else
    % get params matrix
    params = Utilities.getDNActivations(pathStimuli, screeningData.imageIDs, 50, 'fc6');
    options.screenType = 'Object';
    disp('hi');
end

pathOut = [basePath filesep 'PC1PC2Spread' filesep 'NoAxes'];
if ~exist(pathOut)
    mkdir([pathOut]);
end

for cellIndex = 1:length(strctCells)
    psth = screeningData.psth(cellIndex, :);
    timelimits = screeningData.timelimits;
    % note that respLat here has had the timelimit subtracted already
    emptyIndex = cellfun('isempty', screeningData.responses(:, 3));     % Find indices of empty cells
    screeningData.responses(emptyIndex, 3) = {0};                    % Fill empty cells with 0
    cellInQuestion = find(cell2mat(screeningData.responses(:, 3))==strctCells(cellIndex).Name);
    
    if ~isempty(cellInQuestion)
        respLat = screeningData.responses{cellInQuestion, 2};
        cellInfo = strctCells(cellIndex);
        
        prefStimOrder = Utilities.sortByRespMagnitude(screeningData.sortedOrder, screeningData.imageIDs, psth{1, 1}, respLat, screeningData.stimDur);
        topStimOrder = prefStimOrder(1:50);
        bottomStimOrder = prefStimOrder(end-49:end);
        
    end
    
    % make plot
    topStimOrder_params = params(topStimOrder, 1:2);
    bottomStimOrder_params = params(bottomStimOrder, 1:2);
    f = figure;
    hold on
    scatter(topStimOrder_params(:, 1), topStimOrder_params(:, 2), 50, [1 0 0], 'filled');
    scatter(bottomStimOrder_params(:, 1), bottomStimOrder_params(:, 2), 50, [0 0 1], 'filled');
    xline(0);
    yline(0);
    lg = legend({'Top Stimuli', 'Bottom Stimuli'});
    lg.Position = [0.055952385316292,0.194841272111923,0.244642852566072,0.086904759634109];
%     xlim([max([topStimOrder_params(:, 1); bottomStimOrder_params(:, 1)]) -max([topStimOrder_params(:, 1); bottomStimOrder_params(:, 1)])]);
%     ylim([max([topStimOrder_params(:, 2); bottomStimOrder_params(:, 2)]) -max([topStimOrder_params(:, 2); bottomStimOrder_params(:, 2)])]);
    title({[strctCells(cellIndex).brainArea ' ' num2str(strctCells(cellIndex).ChannelNumber) ' ' num2str(strctCells(cellIndex).Name)], 'PC1-PC2 Spread'});
    set(gca, 'visible', 'off');
    set(gca, 'FontWeight', 'bold');
    filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name)];
    print(f, filename, '-dpng', '-r0');
    close all
    
end