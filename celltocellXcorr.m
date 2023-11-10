
%% Script to compute cell-cell xcorrs in different experiments and trial types
% 
% Experiment 1: Screening
%     Trial types - 5 visual categories
%     Compute IT-Amyg, IT-Hipp
%     
% Experiment 2: Imagination
%     Trial types - Screening/encoding, distraction (vis search), Imagination
%     Compute IT-Amyg, IT-Hipp
%     
%     Workflow
%         Grab n IT and m MTL neurons (simultaneously recorded)
%         Per trial type
%             Compute shuffle-corrected xcorr for all nxm cell pairs 
%             Average them
%             
%         Compute xcorr integrals for both sides 

% vwadia Oct2022
%% Setpaths and load in Data

setDiskPaths

Exp = 2; % 1 - screening, 2 = imagination
IT_Amyg = 0;
IT_Hipp = 1;

if Exp == 1 % screening
    
    basePath = [diskPath filesep 'Object_Screening' filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328' filesep '500Stim'];
    taskStruct = load([basePath filesep 'P73CS_ParamObj_Sub_4_Block']); patID = 'P73CS';
    load([basePath filesep 'PsthandResponses'])
    
    % arrange data as long vectors with all collected trial types
    
    
    
elseif Exp == 2 % imagination
    
    basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_3_20210927'];
    taskStruct = load([basePath filesep 'P76CSRec_ReScreen_3_Sub_4_Block']); patID = 'P76CS';
    
    load([basePath filesep 'RecallData_NoFReeRec.mat']) % loads recall data with CRTimeCourse, EncTimeCourse and StrctCells
   
    
end


strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

% arrange data as long vectors with all collected trial types
CROrder = RecallData.CROrder;
[a_CR, b_CR] = sortrows(CROrder);

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';

if IT_Hipp
    
    IT_Hipp = cellfun(@(x) ismember(x, {'LH', 'RH', 'LFFA', 'RFFA'}), strctCELL(:, 4), 'UniformOutput', false);
    strctCells = strctCells(cell2mat(IT_Hipp));
    
    CRTimeCourse = RecallData.CRTimeCourse(cell2mat(IT_Hipp), :);
    
    % sort order - make it easier
    for cellIndex = 1:length(CRTimeCourse)
        
        CRTimeCourse{cellIndex, 1} = CRTimeCourse{cellIndex, 1}(b_CR, :);
        CRTimeCourse{cellIndex, 2} = CRTimeCourse{cellIndex, 2}(b_CR, :);
    end
    
elseif IT_Amyg
    
    IT_Amyg = cellfun(@(x) ismember(x, {'LA', 'RA', 'LFFA', 'RFFA'}), strctCELL(:, 4), 'UniformOutput', false);
    strctCells = strctCells(cell2mat(IT_Amyg));
    
    CRTimeCourse = RecallData.CRTimeCourse(IT_Amyg, :);
    
end



%% compute xcorrs








