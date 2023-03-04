% Toy script to check if Liang's code for a shuffle distribution is being
% computed correctly. 
% 
% Need to run screeningScript first
% vwadia May2021
% 
%%
if strcmp(taskStruct.subID, 'P73CS_Full')
    z_scored = 0;
    score = load('Z:\LabUsers\vwadia\SUAnalysis\ObjectSpace\parameters_2k_synthetic_faces.mat');
    score = score.params;
end

if strcmp(taskStruct.subID, 'P73CS_ParamObj')
    load('Z:\LabUsers\vwadia\SUAnalysis\ObjectSpace\parameters_1593_objects.mat'); % will create score = 1593x50
end
pathOut = [basePath filesep 'STA_and_projections_varWindow'];
if ~exist(pathOut)
    mkdir([pathOut]);
end

for cellIndex = 1:length(strctCells)
    
    if ~isempty(screeningData.avgResponses{cellIndex, 2})
        options.ind_train = screeningData.imageIDs; % use all objects to calculate STA
        options.fam = 0;
        options.unfam = 0;
%         [hfig] = STA_figure_original(screeningData.avgResponses{cellIndex}, score, options); % pass score to this instead of projectedResponses
          [hfig] = face_id_analysis_STA_varun(screeningData.avgResponses{cellIndex}, options); % pass score to this instead of projectedResponses
        
       
        suptitle({['Cell number ' num2str(strctCells(cellIndex).Name)], strctCells(cellIndex).brainArea});
        print(hfig, [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).Name)], '-dpng', '-r0')
        
        
       
        close all;
    end
end