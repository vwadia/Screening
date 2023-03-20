dbstop if error
setDiskPaths

ITCellsOnly = 1;
imageIDs = [1:500]';

%% Combining and recomputing all cells together (takes a while)

% all cells - morning
load([diskPath filesep 'Object_Screening' filesep 'ALLITCells_500Stim_Scrn.mat']);

m_strctCELL = struct2cell(strctCells');
m_strctCELL = m_strctCELL';


m_psths = psths;
m_responses = responses;
m_strctCells = strctCells;

% get rid of non-responsive units 
index = cell2mat(cellfun(@isnan, m_responses(:, 2), 'UniformOutput', false));
m_responses(index(:, 1), :) = [];
m_psths(index(:, 1), :) = [];
m_strctCells(index(:, 1)) = [];

%
% all cells - afternoon. Note that this will reset 'responses'/'psths'/'strctCells'
load([diskPath filesep 'Object_Screening' filesep 'ALLITCells_500Stim_ReScreen.mat']);

a_strctCELL = struct2cell(strctCells');
a_strctCELL = a_strctCELL';


a_psths = psths;
a_responses = responses;
a_strctCells = strctCells;

% get rid of non-responsive units 
index = cell2mat(cellfun(@isnan, a_responses(:, 2), 'UniformOutput', false));
a_responses(index(:, 1), :) = [];
a_psths(index(:, 1), :) = [];
a_strctCells(index(:, 1)) = [];

strctCells = [m_strctCells a_strctCells];
psths = [m_psths; a_psths];
responses = [m_responses(:, 1:3); a_responses(:, 1:3)];

%% adding time field 
% 
% for cm = 1:length(m_strctCells)
%     
%     m_strctCells(cm).time = 'morn';
%     
% end
% 
% for ca = 1:length(a_strctCells)
%     
%     a_strctCells(ca).time = 'aft';
%     
% end

%% computing

load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50


options.screenType = 'Object';
options.ind_train = imageIDs; % use all objects to calculate STA
        
parfor cellIndex = l(strctCells)
    tic
    n_reps = 100;
    p_dist = nan(1, n_reps);
    for dist = 1:n_reps
        [p, pRat] = Utilities.ObjectSpace.linearity_measure_STA(responses{cellIndex, 1}, params, options);
        p_dist(dist) = p;
    end
    
    % save the whole distribution
    strctCells(cellIndex).pvalRamp = p_dist;
    disp(['Finished for cell ' num2str(cellIndex)])
    toc
end

keyboard
% save file
save([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'], 'strctCells', 'psths', 'responses', '-v7.3')

%% adding a session at a time (much faster than re-computing dist for all cells)

% % set new sess path
% newSess = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030'];
% sessTime = 'aft'; % set this manually
% % keyboard
% 
% % load new sess stuff
% n_S = load([newSess filesep 'strctCells']);
% n_PnR = load([newSess filesep 'PsthandResponses']);
% 
% % get rid of non-responsive units
% index = cellfun(@isempty, n_PnR.responses);
% n_PnR.responses(index(:, 1), :) = [];
% n_PnR.screeningPsth(index(:, 1), :) = [];
% n_S.strctCells(index(:, 1)) = [];
% 
% % load params
% load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50
% options.screenType = 'Object';
% options.ind_train = imageIDs; % use all objects to calculate STA
% 
% % compute for new session
% for cellIndex = l(n_S.strctCells)
%     tic
%     n_reps = 100;
%     p_dist = nan(1, n_reps);
%     for dist = 1:n_reps
%         [p, pRat] = Utilities.ObjectSpace.linearity_measure_STA(n_PnR.responses{cellIndex, 1}, params, options);
%         p_dist(dist) = p;
%     end
%     
%     % save the whole distribution
%     n_S.strctCells(cellIndex).pvalRamp = p_dist;
%     n_S.strctCells(cellIndex).time = sessTime;
%     disp(['Finished for cell ' num2str(cellIndex)])
%     toc
% end
% 
% 
% % load all other cells 
% 
% for set = 1:2
%     
%     if set == 1
%         cellSet = 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat';
%     elseif set == 2
%         cellSet = 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat';
%     end
%     % load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])
%     % load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'])
%     load([diskPath filesep 'Object_Screening' filesep cellSet])
%     
%     
%     if strcmp(sessTime, 'morn')
%         % keep all morning cells together
%         times = zeros(length(strctCells), 1);
%         for i = 1:length(strctCells)
%             if strcmp(strctCells(i).time, 'morn')
%                 
%                 times(i, 1) = 1;
%                 
%             end
%             
%         end
%         
%         m_strctCells = [strctCells(times == 1) n_S.strctCells];
%         m_responses = [responses(times == 1, 1:3); n_PnR.responses];
%         m_psths = [psths(times == 1, :); n_PnR.screeningPsth];
%         
%         a_strctCells = strctCells(times == 0);
%         a_responses = responses(times == 0, 1:3);
%         a_psths = psths(times == 0, :);
%         
%         strctCells = [m_strctCells a_strctCells];
%         psths = [m_psths; a_psths];
%         responses = [m_responses(:, 1:3); a_responses(:, 1:3)];
%         
%     elseif strcmp(sessTime, 'aft')
%         % then ust tack onto the end
%         
%         strctCells = [strctCells n_S.strctCells];
%         psths = [psths; n_PnR.screeningPsth];
%         responses = [responses(:, 1:3); n_PnR.responses(:, 1:3)];
%         
%         
%     end
%     
%     % now save new set
%     % save([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'], 'strctCells', 'psths', 'responses', '-v7.3')
%     % save([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'], 'strctCells', 'psths', 'responses', '-v7.3')
%     save([diskPath filesep 'Object_Screening' filesep cellSet], 'strctCells', 'psths', 'responses', '-v7.3')
% 
% end 
% 

%% after manually validating cellmatching per session  - make a new struct

setDiskPaths

dbstop if error
[~, host] = system('hostname');
if strcmp(host(1:end-1), 'DWA644201')
    atCedars = 1;
    xlsFile = 'D:\Dropbox\Caltech\Thesis\Human_work\Cedars\for_Advisors_VarunThesis\IT_ImaginationPaper\Summary of Sessions and Cells Screening+Recall.xlsx';
elseif strcmp(host(1:end-1), 'DESKTOP-LJHLIED')
    atCedars = 0;
    xlsFile = 'E:\Dropbox\Caltech\Thesis\Human_work\Cedars\for_Advisors_VarunThesis\IT_ImaginationPaper\Summary of Sessions and Cells Screening+Recall.xlsx';
else % mac
    atCedars = 0;
    xlsFile = '/Volumes/Macintosh HD/Users/varunwadia/Dropbox/Caltech/Thesis/Human_work/Cedars/for_Advisors_VarunThesis/IT_ImaginationPaper/Summary of Sessions and Cells Screening+Recall.xlsx';
end

% fill this in manually
sheet = 'grabNewCells';
cell_col = 7;
sess_col = 8; % need this because ofredundancy in names
range = [2:120];
% mornInfo = xlsread(xlsFile, sheet, '', 'basic');
[~, ~, mornInfo] = xlsread(xlsFile, sheet, '', 'basic');

mornCellsToExclude = mornInfo(2:end, cell_col); 
cellnums = cellfun(@(x) ~isnan(x), mornInfo(2:size(mornInfo, 1), cell_col));
mornCellsToExclude = cell2mat(mornCellsToExclude(cellnums, 1));

mornSessionsOfExCells = mornInfo(2:end, sess_col);
sessIds = cellfun(@(x) ~isnan(x), mornInfo(2:size(mornInfo, 1), sess_col), 'UniformOutput', false);
fuckthis = cellfun(@(x) length(x) == 1, sessIds, 'UniformOutput', false);
mornSessionsOfExCells = mornSessionsOfExCells(~cell2mat(fuckthis), 1);
% load in all cells 
load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat']);

dupCells = zeros(length(strctCells), 1);
for cellIndex = l(strctCells)
    for c_I = 1:length(mornCellsToExclude)
        if strctCells(cellIndex).Name == mornCellsToExclude(c_I)...
                && strcmp(strctCells(cellIndex).SessionID, mornSessionsOfExCells{c_I})
            dupCells(cellIndex) = 1;
        end
    end  
end

% all cell names
cellNames = cell2mat(responses(:, 3));

% remove morning cells with afternoon pairs 
new_psths = psths(~dupCells, :);
new_strctCells = strctCells(~dupCells);
new_resp = responses(~dupCells, :);

% combine
strctCells = new_strctCells;
responses = new_resp;
psths = new_psths; 

% now save
save([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'], 'strctCells', 'psths', 'responses', '-v7.3')


%% Finally - Saving merged sigramp cells

% load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])
load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'])

idx = zeros(length(strctCells), 1);


for cellIndex = l(strctCells)
    
    p_info = strctCells(cellIndex).pvalRamp;
    
    
    idx(cellIndex) = sum(p_info < 0.01);
end


strctCells = strctCells(idx == length(p_info));
responses = responses(idx == length(p_info), :);
psths = psths(idx == length(p_info), :);

keyboard
save([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500Stim_Scrn_SigRamp.mat'], 'strctCells', 'responses', 'psths', '-v7.3')



%% for targetting - where did we get good yield

setDiskPaths
load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])
% load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'])

idx = zeros(length(strctCells), 1);


for cellIndex = l(strctCells)
    
    p_info = strctCells(cellIndex).pvalRamp;
    
    idx(cellIndex) = sum(p_info < 0.01);
end

s_c = strctCells(idx == 500);

% amke cell arrays 
strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

s_C = struct2cell(s_c');
s_C = s_C';

PT = {};
pt_ids = {'P71CS', 'P73CS', 'P75CS', 'P76CS', 'P77CS_L', 'P77CS_R', 'P78CS', 'P79CS_L', 'P79CS_R', 'P80CS_L', 'P80CS_R'};
pt_ids = fliplr(pt_ids);

for i = 1:length(pt_ids)
    
  PT{i, 1} = pt_ids{i};
  
  if strcmp(pt_ids{i}(end) , 'L')
      
      t1 = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), strctCELL(:, 7)); t2 = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));
      PT{i, 2} = sum(t1+t2 == 2);
      
      t1_sc = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), s_C(:, 7)); t2_sc = cellfun(@(x) strcmp(x, 'LFFA'), s_C(:, 4));
      PT{i, 3} = sum(t1_sc+t2_sc == 2);
      
  elseif strcmp(pt_ids{i}(end) , 'R')
      
      t1 = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), strctCELL(:, 7)); t2 = cellfun(@(x) strcmp(x, 'RFFA'), strctCELL(:, 4));
      PT{i, 2} = sum(t1+t2 == 2); 
      
      t1_sc = cellfun(@(x) strcmp(x, pt_ids{i}(1:end-2)), s_C(:, 7)); t2_sc = cellfun(@(x) strcmp(x, 'RFFA'), s_C(:, 4));
      PT{i, 3} = sum(t1_sc+t2_sc == 2);
      
      
  else
      
      PT{i, 2} = sum(cellfun(@(x) strcmp(x, pt_ids{i}), strctCELL(:, 7)));
      PT{i, 3} = sum(cellfun(@(x) strcmp(x, pt_ids{i}), s_C(:, 7)));    
      
  end
  
  
end  
    
    
