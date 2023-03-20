

%% For imagination related things - see scrn_plot_ramp_summary
%%

% This script is used to combine strctCells and responses for many (screening) sessions AFTER
% the psths, responses, order etc. have been saved.
%% INSTRUCTIONS

% Run screeningScript to get infor per session (respLat etc.) - save those
% Run RecallScript and then plotPerCellRasterandBarPerStim to get strctResp - save those
% load them in here, combine and save

%% for everything else responsive units/sigramp units 

% use loadCellsComputeRampSig 


%% Repeat this for each session  

% % set paths and create strctCells from screeningScript
% % then run this an manually compile

dbstop if error
setDiskPaths

task = 'Recall_Task';
% task = 'ReScreen_Recall'; % collect the screening cells from the recall sessions
% task = 'Object_Screening';

taskPath = [diskPath filesep task];

ITCellsOnly = 1;

%% combine

if strcmp(task, 'Recall_Task') 

    sessID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
        ['Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925'],... % moved ITResponses into screening folder for sess 2
        ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
        ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
        ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731']};    
    
    for s = 1:length(sessID)
        basePath = [diskPath filesep sessID{s}];
        load([basePath filesep 'ITResponses']); % creates strctResp
        s_CRresp{s, 1} = strctResp;
    end
    
elseif strcmp(task, 'ReScreen_Recall') % all afternoon neurons 
    
     sessID = {['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'],...
        ['Recall_Task' filesep 'P76CS' filesep 'RecallScreening_Session_2_20210925'],... % moved ITResponses into screening folder for sess 2
        ['Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_1_20220330'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_2_20220403'],...
        ['Recall_Task' filesep 'P79CS' filesep 'ReScreenRecall_Session_3_20220405'],...
        ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_1_20220728'],...
        ['Recall_Task' filesep 'P80CS' filesep 'ReScreenRecall_Session_2_20220731'],...
        ['Object_Screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030'],...
        ['Object_Screening' filesep 'P82CS' filesep 'ClosedLoopReScreen_Session_1_20230115']};
    
    
    
elseif strcmp(task, 'Object_Screening') % morning sessions only
    sessID = {['Object_Screening' filesep 'P71CS' filesep 'FastObjectScreening_Session_1_20201125'],...%d
        ['Object_Screening' filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328' filesep '500stim'],...%d
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
        ['Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030'],...
        ['Object_Screening' filesep 'P82CS' filesep 'FingerprintScreening_Session_1_20230111'],...
        ['Object_Screening' filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115']}; %d
    
    
end

for s = 1:length(sessID)
        basePath = [diskPath filesep sessID{s}];
        load([basePath filesep 'strctCells']);
        load([basePath filesep 'PsthandResponses']);
        s_cells{s, 1} = strctCells;
        s_resp{s, 1} = responses;
        s_resp{s, 2} = psths;
end


%% once compiled - automatically combine 

% Loads in all the relevant cells across patients and sessions - combines them into a single strctCells
strctCells = Utilities.combineStrcts(s_cells, task);
 
[responses, psths] = Screening.combineResponsesAndPsths(s_resp);
 
if strcmp(task, 'Recall_Task')
    strctResp = Utilities.combineStrcts(s_CRresp, task);
end


%% save cells screening (morning/afternoon sessions)

% IT cells
if ITCellsOnly
    
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    IT_Cells = cellfun(@(x) strcmp(x, 'RFFA') || strcmp(x, 'LFFA'), strctCELL(:, 4));
    
    psths = psths(IT_Cells, :);
    responses = responses(IT_Cells, :);
    strctCells = strctCells(IT_Cells);
    
    if strcmp(task, 'Object_Screening')
        save([taskPath filesep 'AllITCells_500stim_Scrn'], 'responses', 'psths', 'strctCells', '-v7.3');
    elseif strcmp(task, 'ReScreen_Recall')
        save([diskPath filesep 'Object_Screening' filesep 'AllITCells_500stim_ReScreen'], 'responses', 'psths', 'strctCells', '-v7.3');
    elseif strcmp(task, 'Recall_Task')
        save([taskPath filesep 'AllITCells_500stim_Im'], 'responses', 'psths', 'strctCells', '-v7.3'); 
    end
    
else
    
    if strcmp(task, 'Object_Screening')
        %     All possible cells - for autocorr computation
        save([taskPath filesep 'AllCells_500stim_Scrn'], 'responses', 'psths', 'strctCells', '-v7.3');
    elseif strcmp(task, 'ReScreen_Recall')
        save([taskPath filesep 'AllCells_500stim_ReScreen'], 'responses', 'psths', 'strctCells', '-v7.3');  
    elseif strcmp(task, 'Recall_Task')
        save([taskPath filesep 'AllCells_500stim_Im'], 'responses', 'psths', 'strctCells', '-v7.3'); % this is essentially the same as the RescreenRecall cells - except for that one annoying time I didn't get a rescreen and the synth image screen
    end
    
end

%% Go to loadCellsComputeRampSig for sigramp computation (pvalDistribution)


%% Saving cells recall - if you already have all cells
% NOTE: You need to run the first couple of cells of this to get strctResp
setDiskPaths
load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])

idx = zeros(length(strctCells), 1);

for cellIndex = l(strctCells)
    
    if strcmp(strctCells(cellIndex).time, 'aft') && ~strcmp(strctCells(cellIndex).SessionID, 'P81_synth') && ~strcmp(strctCells(cellIndex).SessionID, 'P82CS_CLReScreen')
        idx(cellIndex) = 1;   
    end
    
end

% remove morning (screening) neurons and the CLScreen neurons
strctCells(~idx) = [];
responses(~idx, :) = [];
psths(~idx, :) = [];

%%
% paranoia check - make sure cells are set up correctly
% P76 Session 2 won't be but everything else should be fine
for cellIndex = l(strctCells)
    if ~isequal(strctResp(cellIndex).Name, strctCells(cellIndex).Name)
        keyboard
    end  
end

%%
idx = [];

% now find sigramp neurons 
for cellIndex = l(strctCells)
    
    p_info = strctCells(cellIndex).pvalRamp;  
%     p_info = p_info(1);
    
    idx(cellIndex) = sum(p_info < 0.01);
end

strctCells(idx ~= length(p_info)) = []; 
responses(idx ~= length(p_info), :) = [];
psths(idx ~= length(p_info), :) = [];
strctResp(idx ~= length(p_info)) = [];

save([taskPath filesep 'AllITCells_500Stim_Im_SigRamp'], 'strctCells', 'responses', 'psths', '-v7.3')
save([taskPath filesep 'AllITResponses_500Stim_Im_SigRamp'], 'strctResp', '-v7.3')


%% Old way before merges

% % get rid of non-responsive units 
% index = cellfun(@isempty, responses);
% responses(index(:, 1), :) = [];
% psths(index(:, 1), :) = [];
% strctCells(index(:, 1)) = [];
% 
% if strcmp(task, 'Recall_Task')
%     strctResp(index(:, 1)) = [];
% end
% 
% for cellIndex = l(strctCells)
%     
%     if ~isequal(strctResp(cellIndex).Name, strctCells(cellIndex).Name)
%         keyboard
%     end
%     
%     
% end 
% 
% if strcmp(task, 'Recall_Task')
%     
%     save([taskPath filesep 'AllRespITCells_500stim_Im'], 'responses', 'psths', 'strctCells', '-v7.3');
%     save([taskPath filesep 'AllResponses_RespITCells_500Stim_Im'], 'strctResp', '-v7.3');
%     
% elseif strcmp(task,  'Object_Screening')
%     
%     save([taskPath filesep 'AllRespCells_500stim_Scrn'], 'responses', 'psths', 'strctCells', '-v7.3');
%     
%     
% elseif strcmp(task, 'ReScreen_Recall')
%     
%     save([taskPath filesep 'AllRespCells_500stim_ReScreen'], 'responses', 'psths', 'strctCells', '-v7.3');
%     
% end
