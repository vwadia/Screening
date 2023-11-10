
%% 1. get alexnet activations
lexnet = alexnet; 
% % analyzeNetwork(net);
% % for activations function
threeD = 1;
grayImages = Utilities.getImageDescriptions(pathStimuli, 227, screeningData.imageIDs, threeD);
% %%
% load([basePath filesep '3DImages'])
% make this a function
act = activations(lexnet, grayImages, 'fc6', 'OutputAs', 'rows'); 
% r_1593 = act;

respScreeningImagesLayer = act;
% respScreeningImages = respScreeningImagesLayer;

% ctr = 1;
% respScreeningImagesLayer = respScreeningImages{ctr};
[coeff,score,latent,tsquared,explained,mu] = pca(respScreeningImagesLayer);
score = score(:, 1:50); % you can do this in STA figure or here
params = score;

%% Same as cell 1 ^

% my explained, liangs explained, explanable
screeningData.explainedVar = zeros(length(strctCells), 2);
ctr = 1;
if strcmp(taskStruct.subID, 'P73CS_ParamObj')
    screenType = 'Object';    
    load([diskPath filesep 'ObjectSpace\parameters_1593_objects.mat']); % will create params = 1593x50
elseif strcmp(taskStruct.subID, 'P73CS_Full')
    screenType = 'Face';
    load([diskPath filesep 'ObjectSpace\parameters_2k_synthetic_faces.mat']); % will create shape-appearance parameters. params = 2000x50
    params = params(1:667, :);
end

%% for FR prediction per cell - validating axis model
for cellIndex = l(strctCells)   
    if ~isempty(screeningData.responses{cellIndex, 2}) % responsive cells only
        
        % observed responses
        resp = screeningData.responses{cellIndex};
        
        % predicted responses
        [pred_resp, obs_resp] = Utilities.computePredictedResponses(resp, params, screeningData.imageIDs, screenType);
        screeningData.pred_resp{cellIndex} = pred_resp;
        
        % calculate explained variance
        screeningData.explainedVar(cellIndex, 1) = Utilities.computeExplainedVariance(obs_resp, pred_resp);
        
    end    
end

%% Liangs code for explained variance per cell
% 
% for cellIndex = l(strctCells)   
%     % responsive cells only
%     if ~isempty(screeningData.responses{cellIndex, 2}) 
%         
%         % collect responses
%         resp = screeningData.responses{cellIndex};
%         % calculate ev
%         [ev] = STA_sub_cross_val(resp, params, 'sta', 0.01); 
%         screeningData.explainedVar(cellIndex, 2) = ev;
%     end 
% end


%% Explanable variance

screeningData.explainableVar = zeros(length(strctCells), 1);

for cellIndex = l(strctCells)
    
    if ~isempty(screeningData.responses{cellIndex, 2}) % responsive cells only
        respLat = screeningData.responses{cellIndex, 2};
        screeningData.explainableVar(cellIndex, 1) = Utilities.computeExplainableVariance(screeningData.psth(cellIndex, :),...
            screeningData.sortedOrder, respLat, screeningData.stimDur, -screeningData.timelimits(1)*1e3);              
    end
    
end

%% comparison figure

f = figure; 
hold on; 
scatter(screeningData.explainedVar(:, 1), screeningData.explainableVar); 
xlabel('Explained variance 50d space')
ylabel('Explainable variance (split half)')
title('1593 Object Screen P73CS FFA')
% title('667 Face Screen P73CS FFA')
% title('AIC Screen P73CS FFA')

% xlim([-0.035 0.035]);
xline(0);
yline(0);
% print(f, [basePath filesep 'explained_explanable_variance_scatter'],'-dpng', '-r0' ) 

%% compute variance per image - checking regression model fit for image decoding 

% where tf does pred_feat come from?
for image = 1:size(pred_feat, 1)
    
        % observed responses
        obs_resp = params(image, :);
        % predicted responses
        pred_resp = pred_feat(image, :);
        
        % calculate explained variance
        screeningData.eVImage(image, 1) = computeExplainedVariance(obs_resp, pred_resp);

end





                


