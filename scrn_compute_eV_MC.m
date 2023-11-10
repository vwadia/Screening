
%% Comparing the explained variance of IT cells responses to both faces and objects - across multiple models

%% workFlow

% Load in cells and responses
% load in stimuli in a cell array (For Alexnet, do a redundancy check by selecting out of parameters_1593 and creating new descriptions = assert equality)
% Load in pretrained networks and get unit activations for the images, do pca and create compact representation
% Use LOO/10-fold cross validation to get a predicted response to each image per neuron
% 1. compute eV and average across all cells 2. compute percentage of explanable variance found 3. ev divded by noise ceiling

% compute encoding error (angle b/w pred resp vec per image (cellsX1) and target/distractor)
% compute decoding error (Same except pred feat vec per image (numFeatX1)) - is this correct?

% Note: Also need to redo procedure for all layers of each NN

setDiskPaths
taskPath = 'Object_Screening';
addpath(genpath('ObjectSpace'))

% For Objects

screenType = 'Object';

% 64 cells, 500 objects
% Networks: Alexnet, VGG-19, VGG-16, VGG-Face, CORnet(4 versions), 2D morph (faces)

% load in the neural data
% Responsive cells
% load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat']); RespCells = true; SigCells = false; layermat = 'fc7'; layerpy = 'IT_output'; stimDir = '500Stimuli'; imageIDs = [1:500]';


% load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50
% options.screenType = screenType;
% for cellIndex = 1:length(strctCells)
%     
%     options.ind_train = imageIDs; % use all objects to calculate STA
%     [p, pRat] = Utilities.ObjectSpace.linearity_measure_STA(responses{cellIndex, 1}, params, options);
%     strctCells(cellIndex).pvalRamp = p;
% 
% end
% strctCELLS = reshape(struct2cell(strctCells), [length(fieldnames(strctCells)) length(strctCells)]);
% fields = fieldnames(strctCells);
% sig = cellfun(@(x) strcmp('pvalRamp', x), fields);
% pvals_ramp = cell2mat(strctCELLS(sig, :)');
% clearvars params


% sig ramp cells
load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500Stim_Scrn_SigRamp.mat']); SigCells = true; RespCells = false; layermat = 'fc7'; layerpy = 'IT_output'; stimDir = '500Stimuli'; imageIDs = [1:500]';
% load([diskPath filesep 'Object_Screening' filesep 'old_cellStructs' filesep 'MergedITCells_500Stim_Scrn_SigRamp.mat']); SigCells = 1; layermat = 'fc7'; layerpy = 'IT_output'; stimDir = '500Stimuli'; imageIDs = [1:500]';
layermat = 'comparison'; layerpy = 'comparison_output'; 

% load([diskPath filesep taskPath filesep 'AllFaceCells_1593_Objects.mat']); FaceCells = 1; layermat = 'fc7'; layerpy = 'IT_output'; stimDir = '201Faces'; imageIDs = [1:201]';
siphonOffFace = false;

% siphon off face responses (if  face cells)
if siphonOffFace
    faceOrder = sortedOrder(catOrder == 1);
    faceInds = unique(faceOrder);
    responses(:, 1) = cellfun(@(x) x(faceInds, 1), responses(:, 1), 'UniformOutput', false);
    psths(:, 1:2) = cellfun(@(x) x(catOrder == 1, :), psths(:, 1:2), 'UniformOutput', false);
    strctCells = faceCells;
    sortedOrder = faceOrder-411;
end

if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
    numMdls = 8;
    numLayers = 4;
else
    numMdls = 9;
    numLayers = 1;
end

explainedVar = zeros(length(strctCells), numMdls); % neurons x models
explainableVar = zeros(length(strctCells), 1);

divideByNoiseCeiling = false; % don't compute explainable variance if true

%%
ndim = 50; % number of features to use
for modelNum = 1:numMdls
    
    % load in model features
    switch modelNum
        case 1
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
        case 2
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_vggface_' layermat '_' stimDir '.mat']]);
        case 3
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_vgg16_' layermat '_' stimDir '.mat']]);
        case 4
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_vgg19_' layermat '_' stimDir '.mat']]);
        case 5
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_CORnetZ_' layerpy '_' stimDir '.mat']]);
        case 6
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_CORnetS_' layerpy '_' stimDir '.mat']]);
        case 7
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_CORnetRT_' layerpy '_' stimDir '.mat']]);
        case 8
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_CORnetR_' layerpy '_' stimDir '.mat']]);
        case 9
            load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_eigen_' stimDir '.mat']]);
            
    end
     if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
         allfeats = params;
     end
    
    % get predicted response for each image per neuron
%     for cellIndex = 1:length(strctCells)
%           
%         observed responses
%         resp = responses{cellIndex, 1};
%         
%         predicted responses
%         [pred_resp{cellIndex}, obs_resp] = Utilities.computePredictedResponses(resp, params, imageIDs, screenType);
%         
%         calculate explained variance
%         explainedVar(cellIndex, modelNum) = Utilities.computeExplainedVariance(obs_resp, pred_resp{cellIndex});
%         
%     end
   
    
    for nl = 1:numLayers
        if numLayers > 1
            params = allfeats{nl};
        end
        params = params(:, 1:ndim);
        for cellIndex = 1:length(strctCells)
            
            % collect responses
            resp = responses{cellIndex, 1};
            
            % calculate ev
%             method = 'linear_regression';
            method = 'sta';
            [ev] = Utilities.ObjectSpace.STA_sub_cross_val(resp, params, method, 0.01);
            explainedVar(cellIndex, modelNum, nl) = ev;
            
        end
    end
    
    
    
end
explainedVar_copy = explainedVar;

%% noise ceiling 

A = cellfun(@(x) x', psths(:, 1), 'UniformOutput', false);


for cellIndex = 1:length(strctCells)
    ras = A{cellIndex, 1};
    if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
        order = repelem(imageIDs, 6);
    elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
        order3 = load([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917' filesep 'PsthandResponses']); % had a few trials cut off
        order = order3.order;
    else
        order = repelem(imageIDs, 4);
    end
    for im = 1:length(imageIDs)
        stimRas = ras(:, find(order == imageIDs(im)));
        A{cellIndex, im} = stimRas;
    end    
end
% ctr = 1;
for i=1:size(A,1)
    for j=1:size(A,2)
      hist=A{i,j}; 
      if strcmp(strctCells(i).SessionID, 'P71CS_Fast')
          stimDur = 167;
      else
          stimDur = 267;
      end
      if responses{i, 2}+stimDur > size(hist, 1)
          B{i,j}=mean(hist(responses{i, 2}:end,:),1);% Choose a time window you want to accumulate the spikes into a firing rate 
%           disp(i);
%           ctr = ctr + 1;
      else
          B{i,j}=mean(hist(floor(responses{i, 2}):floor(responses{i, 2})+stimDur,:),1);% Choose a time window you want to accumulate the spikes into a firing rate
      end
      fir(i,j)=mean(B{i,j});  %mean firing rate
    end
end
%compute noise ceiling
for i=1:size(B,1)
    for j=1:size(B,2)
        rast=B{i,j}';
        num(i,j)=length(rast);
        if length(rast)>1
            meu(i,j)=std(rast,0,1)/sqrt(length(rast));
        end
    end
    ev00(i)=1-sum(meu(i,num(i,:)>1).^2)/sum((fir(i,num(i,:)>1)-mean(fir(i,num(i,:)>1))).^2);
end
evn=ev00;



%% explainable variance 
% does it make sense for this to ever be negative? maybe if the responses
% are just noise

if ~divideByNoiseCeiling
    if strcmp(stimDir, '500Stimuli') && SigCells
        load([diskPath filesep taskPath filesep 'ExV_500Stim_1000Reps_SigRampCells']);
    elseif strcmp(stimDir, '500Stimuli') && RespCells
        keyboard % haven't recomputed this for 218 cells
%         load([diskPath filesep taskPath filesep 'ExV_500Stim_1000Reps_RespCells']);
    else
%       % n_reps = 1000 takes 8hrs on laptop for 235 cells, 2 hr 15 min for on 4 cores for 218
        % 
        tic
        n_reps = 1000;
        n_cells = length(strctCells);
        parfor rep = 1:n_reps
            for cellIndex = 1:n_cells%length(strctCells)
                % calculate explainable variance
                respLat = round(responses{cellIndex, 2});
                if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
                    sortedOrder = repelem(imageIDs, 6);
                    stimDur = 166.6250;
                    timelimits = [-0.17, 0.33];
                elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
                    order3 = load([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917' filesep 'PsthandResponses']); % had a few trials cut off
                    sortedOrder = order3.order;
                    stimDur = 266.6250;
                    timelimits = [-0.17, 0.53];
                else
                    sortedOrder = repelem(imageIDs, 4);
                    %         sortedOrder = faceOrder-411;
                    stimDur = 266.6250;
                    timelimits = [-0.17, 0.53];
                end
                explainableVar(cellIndex, rep) = Utilities.computeExplainableVariance(psths(cellIndex, :),...
                    sortedOrder, respLat, stimDur, -timelimits(1)*1e3);
            end
        end
        toc
    end
    expble_var = mean(explainableVar, 2);
    
    % if exist('RespCells', 'var') && RespCells == 1
    %     save([diskPath filesep taskPath filesep 'ExV_500Stim_1000Reps_RespCells'], 'explainableVar', '-v7.3')
    %
    % elseif exist('SigCells', 'var') && SigCells == 1
%         save([diskPath filesep taskPath filesepr 'ExV_500Stim_1000Reps_SigRampCells'], 'explainableVar', '-v7.3')
    %
    % end
    
    
    % save([diskPath filesep taskPath filesep 'ExV_500Stim_1000Reps_RespCells'], 'explainableVar', '-v7.3')
    % save([diskPath filesep taskPath filesep 'ExV_500Stim_1000Reps_SigRampCells'], 'explainableVar', '-v7.3') 
else
    
    % make xpble_var the noise ceiling
    % This is stevens way
    expble_var = evn';
    
end


%% 

if divideByNoiseCeiling
    highnc = (evn > 0.1)';
    tokeep = highnc;
else %use both expV and NC
    highev = expble_var >= 0.1;
    %     highnc = (evn > 0.1)';
    %     tokeep = (highnc + highev)  == 2;
    tokeep = highev;

end

expble_var = expble_var(tokeep);
explainedVar = explainedVar(tokeep, :, :);

keptCells = strctCells(tokeep);

% subsample - test
% sub_sample = randi(sum(tokeep), [floor(sum(tokeep)/2) 1]);
% explainedVar = explainedVar(sub_sample, :);
% expble_var = expble_var(sub_sample, 1);


%% 
% ratio for ev/expv for cells with high noise ceiling
% adj_eV = explainedVar(evn > 0.1, :)./expble_var(evn > 0.1, :);

adj_eV = explainedVar./expble_var;
adj_eV(adj_eV > 1) = 1;
% adj_eV(adj_eV < -1) = -0.5;
if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
    err = std(adj_eV)./sqrt(size(adj_eV, 1));
    err = reshape(err, [numMdls numLayers]); % figure this out for grouped plots
%     err2 = [];
%     for i = 1:numMdls
%         err2 = [err2 err(i, :)];
%     end
%     err = err2;
    adj_eV =  reshape(mean(adj_eV, 1), [numMdls numLayers]);
    
    X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R'});
    X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R'});
else
    
    X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    
    
    err = std(adj_eV)./sqrt(size(adj_eV, 1));
    
end

% adj_eV = [adj_eV(:, 1:6) adj_eV(:, 8)];
% X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'eigenmodel'});
% X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'eigenmodel'});
%%
f = figure;

hold on

if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
    set(gcf,'Position',get(0,'Screensize'))

    b = bar(X, adj_eV);
    set(gca, 'FontSize', 20, 'FontWeight', 'bold');

%     sgtitle({['Ev/Expv across layers of models'], [num2str(ndim) ' features, obj responses - sigramp cells (sta)']});
if divideByNoiseCeiling
    sgtitle({['Ev/NC across layers of models'], [num2str(ndim) ' features']}, 'FontSize', 20, 'FontWeight', 'bold');
    filename = [diskPath filesep taskPath filesep 'forPaper' filesep 'eV_divNC' num2str(ndim) 'featsLayerComparison_Objects_' method '_forPaper'];
    
else
    sgtitle({['Ev/Expv across layers of models'], [num2str(ndim) ' features']}, 'FontSize', 20, 'FontWeight', 'bold');
    filename = [diskPath filesep taskPath filesep 'forPaper' filesep 'eV_' num2str(ndim) 'featsLayerComparison_Objects_' method '_forPaper'];
    
end
%     sgtitle({['Explained variance across layers of models'], [num2str(ndim) ' features']}, 'FontSize', 20, 'FontWeight', 'bold');
    
    % putting error bars correctly
    [ng, nb] = size(adj_eV);
    endPts = nan(nb, ng);
    for j = 1:nb
        endPts(j, :) = b(j).XEndPoints;
    end
    errorbar(endPts', adj_eV, err, 'k', 'linestyle', 'none')
else
    b = bar(X, mean(adj_eV, 1));
%     b.xticklabels = ''
    errorbar(X, mean(adj_eV, 1), err, 'k', 'LineStyle', 'none')
    
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
%     sgtitle({['Ratio of explained to explainable variance in neurons with high NC'], [num2str(ndim) ' features, obj responses - sigramp cells (sta)']}, 'FontWeight', 'bold');
if divideByNoiseCeiling
    filename = [diskPath filesep taskPath filesep 'forPaper' filesep 'eV_divNC_' num2str(ndim) 'featsModelComparison_Objects_' method '_forPaper'];
    
    sgtitle({['Ev/NC across models'], [num2str(ndim) ' features']}, 'FontWeight', 'bold');
else
    filename = [diskPath filesep taskPath filesep 'forPaper' filesep 'eV_' num2str(ndim) 'featsModelComparison_Objects_' method '_forPaper'];
    
    sgtitle({['Ev/Expv across models'], [num2str(ndim) ' features']}, 'FontWeight', 'bold');
end
%     sgtitle({['Explained variance across models'], [num2str(ndim) ' features']}, 'FontWeight', 'bold');
    
end

if exist('RespCells', 'var') && RespCells == 1
    filename = [filename '_RespCells'];
elseif exist('SigCells', 'var') && SigCells == 1
    filename = [filename '_SigRamp'];
end

% filename = [diskPath filesep taskPath filesep 'eV_comparison_forHSNPoster']
% text(X(2), 0.375, 'fc6', 'color', [0 0.4470 0.7410], 'FontSize', 18, 'FontWeight', 'bold' )
% text(X(2), 0.36, 'relu6', 'color', [0.8500 0.3250 0.0980], 'FontSize', 18, 'FontWeight', 'bold' )
% text(X(2), 0.345, 'fc7', 'color', [0.9290 0.6940 0.1250], 'FontSize', 18, 'FontWeight', 'bold' )
% text(X(2), 0.33, 'relu7', 'color', [0.4940 0.1840 0.5560], 'FontSize', 18, 'FontWeight', 'bold' )
% 
% text(X(7), 0.375, 'V1', 'color', [0 0.4470 0.7410], 'FontSize', 18, 'FontWeight', 'bold' )
% text(X(7), 0.36, 'V2', 'color', [0.8500 0.3250 0.0980], 'FontSize', 18, 'FontWeight', 'bold' )
% text(X(7), 0.345, 'V4', 'color', [0.9290 0.6940 0.1250], 'FontSize', 18, 'FontWeight', 'bold' )
% text(X(7), 0.33, 'IT', 'color', [0.4940 0.1840 0.5560], 'FontSize', 18, 'FontWeight', 'bold' )
% print(f, filename, '-dpng', '-r0')

%% Histogram
%% Use the responsive cells - not the sig ramp cells

idx = zeros(length(keptCells), 1);
for cellIndex = 1:length(keptCells)
    
    p_info = keptCells(cellIndex).pvalRamp;
      
    idx(cellIndex) = sum(p_info < 0.01);
end
pvals_ramp = idx ~= length(p_info);



% pvals_ramp = pvals_ramp(tokeep);
h1 = adj_eV(pvals_ramp < 0.01, 1);
h2 = adj_eV(pvals_ramp >= 0.01, 1);
% histogram of ev/exv
f = figure; 
hold on; 
histogram(h1, -0.2:0.1:1, 'FaceColor', [0.6350 0.0780 0.1840], 'EdgeColor', [0 0 0]); 
histogram(h2, -0.2:0.1:1, 'FaceColor', [0.4940 0.1840 0.5560], 'EdgeColor', [0 0 0]);

title('Explained variance of IT neurons', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Explained variance')
ylabel('Number of cells')
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% ylim([0 8]);
if exist('RespCells', 'var') && RespCells == 1
    lgnd = legend({'Sig ramps', 'Non-sig ramps'}, 'FontWeight', 'bold');

    filename = [diskPath filesep taskPath filesep 'Hist_eV_respCells_AlexNet'];
elseif exist('SigCells', 'var') && SigCells == 1
    lgnd = legend({'Sig ramps'}, 'FontWeight', 'bold');

    filename = [diskPath filesep taskPath filesep 'Hist_eV_SigCells_AlexNet'];
end

% print(f, filename, '-dpng', '-r0');
% for i = 1:4
%     test{i} = explainedVar(:, :, i);
% end
% 
% test_2 = cellfun(@(x) mean(x, 1), test, 'UniformOutput', false);
% print(f, filename, '-dpng', '-r0')

%% For Faces

% 63 cells, 77 faces 
%     AND
% 21 face cells, 201 faces
% Networks: 2D morph, Eigenface, AlexNet, VGG-Face, VGG-19, VGG-16, CORnet(3 versions), Hebbian learning? 

%% Making params from model outputs

