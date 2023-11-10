%% 
%% Column 1 - my analysis, column 2 - Liangs, column 3 - Yuelins
screeningData.explainedVar = zeros(length(strctCells), 2);

if strcmp(taskStruct.subID, 'P73CS_ParamObj')
    if ~exist('params')
        load([diskPath filesep 'ObjectSpace\parameters_1593_objects.mat']); % will create params = 1593x50
    end
    params = score;
elseif strcmp(taskStruct.subID, 'P73CS_Full')
    if ~exist('params')
        load([diskPath filesep 'ObjectSpace\parameters_2k_synthetic_faces.mat']); % will create shape-appearance parameters. params = 2000x50
        params = params(1:667, :);
    end
end

%% Liangs code for explained variance per cell

for cellIndex = l(strctCells)   
    % responsive cells only
    if ~isempty(screeningData.responses{cellIndex, 2}) 
        
        % collect responses
        resp = screeningData.responses{cellIndex};
        % calculate ev
        [ev] = STA_sub_cross_val(resp, params, 'sta', 0.01); 
        screeningData.explainedVar(cellIndex, 1) = ev;
    end 
end


%% Yuelin's code for the same
%% 1. get alexnet activations
lexnet = alexnet; 
% analyzeNetwork(net);
% for activations function
threeD = 1;
grayImages = getImageDescriptions(pathStimuli, 227, screeningData.imageIDs, threeD);

load([basePath filesep '3DImages'])
% make this a function
act = activations(lexnet, grayImages, 'fc6', 'OutputAs', 'rows'); 
respScreeningImagesLayer{1} = act
respScreeningImages = respScreeningImagesLayer;

% % to manually get activations
grayImages_single = getImageDescriptions(pathStimuli, 227, screeningData.imageIDs, 0);
%  
% caffNet = load('imagenet-caffe-alex.mat'); % this has the correct format
% respScreeningImages_2 = getLayerActivations(caffNet, grayImages_single);

% save([basePath filesep '3DImages'], 'grayImages', '-v7.3')




