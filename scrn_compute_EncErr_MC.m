% same recipe as ev calculation
% except you zscore the responses first
% save the predicted responses per image so you have a population vector
% compare the angle of the predicted population response vector
% with that of the target population vector and the pop vector for a distractor face

% if the angle between the pred and target > angle between the pred and distractor = 1, else 0

%% setpaths

setDiskPaths
taskPath = 'Object_Screening';
addpath(genpath('ObjectSpace'))

%% load in data
screenType = 'Object';

load([diskPath filesep taskPath filesep 'MergedITCells_500Stim_Scrn_SigRamp.mat']); layermat = 'fc7'; layerpy = 'IT_output'; stimDir = '500Stimuli'; imageIDs = [1:500]';
% layermat = 'comparison'; layerpy = 'comparison_output';

if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
    numMdls = 8;
    numLayers = 4;
else
    numMdls = 9;
    numLayers = 1;
end


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

for i=1:size(A,1)
    for j=1:size(A,2)
        hist=A{i,j};
        B{i,j}=mean(hist(responses{i, 2}:responses{i, 2}+267,:),1);% Choose a time window you want to accumulate the spikes into a firing rate
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
% responses = responses(evn > 0.1, :);

%%
tic
n_dist = 1; % 51 min for 10 dist 50 dim
n_reps = 1000;
ndim = 50; % number of features to use
enc_acc_overall = zeros(numMdls, numLayers, n_dist);

% ~ 5min across models
% ~ 20min across layers
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
    else
        params = params(:, 1:ndim);
    end
    
   
    for nl = 1:numLayers
        if numLayers > 1
            params = allfeats{nl};
        end
        params = params(:, 1:ndim);
        
         % get predicted response for each image per neuron
        for cellIndex = 1:length(strctCells)
            
            % observed responses - note zscore
            resp = zscore(responses{cellIndex, 1});
            
            % predicted responses
%             method = 'linear_regression';
            method = 'sta';
            [pred_resp(:, cellIndex), obs_resp(:, cellIndex)] = Utilities.computePredictedResponses(resp, params, imageIDs, screenType, method);
            
        end
        
        pred_resp = pred_resp(:, evn > 0.1);
        obs_resp = obs_resp(:, evn > 0.1);
        % project prediction onto target and a randomly chosen distractor
        for nd = 1:n_dist
            enc_acc = zeros(length(imageIDs), n_reps);
            for s = 1:size(pred_resp, 1) % number of images
                
                for rep = 1:n_reps
                    
                    img = pred_resp(s, :); % img
                    target = obs_resp(s, :); % target
                    
                    allInds = 1:size(pred_resp, 1);
                    distInds = setdiff(allInds, s);
                    dist = obs_resp(randsample(distInds, nd), :); % distractor(s)
                    
                    % projection decoding - why are these different???
                    if sum(img*target' > img*dist')/nd == 1
                        %                 enc_acc(s, rep) = sum(img*target' > img*dist')/n_dist;
                        enc_acc(s, rep) = 1;
                    end
                    
                    % nn decoding - why are these different???
%                     distance = [];
%                     for run = 1:nd
%                         distance(run) = norm(dist(run, :) - img);
%                     end
%                     
%                     distance = [distance norm(target - img)];
%                     [min_dist, min_idx] = min(distance);
%                     if min_idx == length(distance)
%                         enc_acc(s, rep) = 1;
%                     end
                    
                end
            end
            enc_acc_per_stim(:, modelNum, nl) = mean(enc_acc, 2);
            enc_acc_overall(modelNum, nl, nd) = mean(mean(enc_acc));
        end
    end
end
toc
%%
if n_dist == 1
    X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    % enc_err = 1-enc_acc_overall(:, 1);
    enc_err = 1-enc_acc_overall;
    
    err = std(enc_acc_per_stim)./sqrt(length(imageIDs));
    if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
        X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R'});
        X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R'});
        
        err = reshape(err, [8 4]);
    end
    
    f = figure;
    hold on
    b = bar(X, enc_err);
    
    % place error bars correctly
    if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
        set(gcf,'Position',get(0,'Screensize'))
        set(gca, 'FontSize', 20, 'FontWeight', 'bold');
        
        [ng, nb] = size(enc_err);
        endPts = nan(nb, ng);
        for j = 1:nb
            endPts(j, :) = b(j).XEndPoints;
        end
        errorbar(endPts', enc_err, err, 'k', 'linestyle', 'none')
        %     filename = [diskPath filesep taskPath filesep 'EncErr_100dimsLayerComparison_Objects_sta_highNC_LOO_VC']
        filename = [diskPath filesep taskPath filesep 'EncErr_' num2str(ndim) 'dimsLayerComparison_Objects_' method '_forPaper'];
        sgtitle(['Encoding error across Layers ' num2str(ndim) ' features'], 'FontSize', 20, 'FontWeight', 'bold');
        
    else
        errorbar(X, enc_err, err, 'k', 'LineStyle', 'none');
        %     filename = [diskPath filesep taskPath filesep 'EncErr_100dimsModelComparison_Objects_sta_highNC_LOO_VC']
        filename = [diskPath filesep taskPath filesep 'EncErr_' num2str(ndim) 'dimsModelComparison_Objects_' method '_forPaper'];
        sgtitle(['Encoding error across models ' num2str(ndim) ' features'], 'FontWeight', 'bold');
        set(gca, 'FontSize', 14, 'FontWeight', 'bold');
        
        
    end
else
    
    
    Fontsize = 20;
    % enc_acc_overall = dco_copy; n_dist = 10;
    X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    cols = Utilities.distinguishable_colors(numMdls);
    f = figure;
    set(gcf,'Position',get(0,'Screensize'))
    hold on
    for j = 1:numMdls
        plot(1-enc_acc_overall(j, :), 'Color', cols(j, :), 'LineWidth', 1.5)
    end
    plot(1-(1./[2:n_dist+1]), '--k', 'LineWidth', 2)
    lgnd = legend(X);
    lgnd.FontSize = 16;
    lgnd.Position = [0.908,0.676,0.089,0.249];
    filename = [diskPath filesep taskPath filesep 'EncErr_' num2str(ndim) 'LinePlot_' method '_highNC_VC'];
    sgtitle(['Encoding error across models '  num2str(ndim) ' features'], 'FontSize',Fontsize, 'FontWeight', 'bold');
    ylabel('Encoding error','FontSize',Fontsize, 'FontWeight', 'bold')
    xlabel('Number of distractors','FontSize',Fontsize, 'FontWeight', 'bold')
    set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
    
end
    print(f, filename, '-dpng', '-r0')

