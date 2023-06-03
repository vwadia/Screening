% Using responses to faces to predict the model features% same recipe as ev calculation
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

fr_mat = [];
% create fr matrix from responses
for cnt = 1:length(responses)
    fr_mat = [fr_mat responses{cnt, 1}];
end



%% noise ceiling

A = cellfun(@(x) x', psths(:, 1), 'UniformOutput', false);

for dimIndex = 1:length(strctCells)
    ras = A{dimIndex, 1};
    if strcmp(strctCells(dimIndex).SessionID, 'P71CS_Fast')
        order = repelem(imageIDs, 6);
    elseif strcmp(strctCells(dimIndex).SessionID, 'P76CSRec_ReScreen')
        order3 = load([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917' filesep 'PsthandResponses']); % had a few trials cut off
        order = order3.order;
    else
        order = repelem(imageIDs, 4);
    end
    for im = 1:length(imageIDs)
        stimRas = ras(:, find(order == imageIDs(im)));
        A{dimIndex, im} = stimRas;
    end
end

for i=1:size(A,1)
    for j=1:size(A,2)
        hist=A{i,j};
        B{i,j}=mean(hist(floor(responses{i, 2}):floor(responses{i, 2})+267,:),1);% Choose a time window you want to accumulate the spikes into a firing rate
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
% fr_mat = fr_mat(:, evn>0.1);

%%
tic
n_dist = 1; % 50 min for 10 dist and 50 dim
n_reps = 1000;
ndim = 50; % number of features to use
dec_acc_overall = zeros(numMdls, numLayers, n_dist);

% ~5 mins across models
% ~20 mins across models
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
    
    
    for nl = 1:numLayers
        if numLayers > 1
            params = allfeats{nl};
        end
        params = params(:, 1:ndim);
        
        % zscoring features
        params = zscore(params);
        
        % get predicted response for each image per neuron
        for dimIndex = 1:ndim
            
            % observed responses - note zscore
            feat = params(:, dimIndex);
            
            % predicted features
            % currently I am not normalizing each img fr vector per cell...is
            % this an issue?
%             method = 'linear_regression';
            method = 'sta';
            [pred_feat(:, dimIndex), obs_feat(:, dimIndex)] = Utilities.computePredictedNNFeatures(feat, fr_mat, imageIDs, method);
            
        end
        
        
        % project prediction onto target and a randomly chosen distractor
        for nd = 1:n_dist
            dec_acc = zeros(length(imageIDs), n_reps);
            for s = 1:size(pred_feat, 1) % number of images
                
                for rep = 1:n_reps
                    
                    img = pred_feat(s, :); % img
                    target = obs_feat(s, :); % target
                    
                    allInds = 1:size(pred_feat, 1);
                    distInds = setdiff(allInds, s);
                    d_i = randsample(distInds, nd);
                    dist = obs_feat(d_i, :); % distractor(s)
                    
                    if sum(img*target' > img*dist')/nd == 1
                        dec_acc(s, rep) = 1;
                    end
                    
                    % nn decoding - why are these different???
                    %                 distance = [];
                    %                 for run = 1:nd
                    %                     distance(run) = norm(dist(run, :) - img);
                    %                 end
                    %
                    %                 distance = [distance norm(target - img)];
                    %                 [min_dist, min_idx] = min(distance);
                    %                 if min_idx == length(distance)
                    %                     dec_acc(s, rep) = 1;
                    %                 end
                    
                end
            end
            dec_acc_per_dim(:, modelNum, nl) = mean(dec_acc, 2);
            dec_acc_overall(modelNum, nl, nd) = mean(mean(dec_acc));
        end
    end
end
toc

%%
if n_dist == 1
    X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    % dec_err = 1-dec_acc_overall(:, 1);
    dec_err = 1-dec_acc_overall;
    
    err = std(dec_acc_per_dim)./sqrt(length(imageIDs));
    if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
        X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R'});
        X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R'});
    
        err = reshape(err, [8 4]);
    end
    
    f = figure;
    hold on
    b = bar(X, dec_err);
    
    % place error bars correctly
    if strcmp(layermat, 'comparison') && strcmp(layerpy, 'comparison_output')
        set(gcf,'Position',get(0,'Screensize'))
        set(gca, 'FontSize', 20, 'FontWeight', 'bold');
        [ng, nb] = size(dec_err);
        endPts = nan(nb, ng);
        for j = 1:nb
            endPts(j, :) = b(j).XEndPoints;
        end
        errorbar(endPts', dec_err, err, 'k', 'linestyle', 'none')
        sgtitle(['Decoding error across Layers ' num2str(ndim) ' features'], 'FontSize', 20, 'FontWeight', 'bold');
        filename = [diskPath filesep taskPath filesep 'DecErr_' num2str(ndim) 'featsLayerComparison_Objects_' method '_forPaper'];
    
    else
        errorbar(X, dec_err, err, 'k', 'LineStyle', 'none');
        sgtitle(['Decoding error across models ' num2str(ndim) ' features'], 'FontWeight', 'bold');
        filename = [diskPath filesep taskPath filesep 'DecErr_' num2str(ndim) 'featsModelComparison_Objects_' method '_forPaper'];
        set(gca, 'FontSize', 14, 'FontWeight', 'bold');

    end
    
else
    % line plot
    Fontsize = 20;
    X = categorical({'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    X = reordercats(X, {'alexnet', 'vggface', 'vgg16', 'vgg19', 'CORNet-Z', 'CORNet-S', 'CORNet-RT', 'CORNet-R', 'eigenmodel'});
    
    % enc_acc_overall = dco_copy; n_dist = 10;
    cols = Utilities.distinguishable_colors(numMdls);
    % dec_err = reshape(dec_err, [numMdls n_dist]);
    dec_acc_overall = reshape(dec_acc_overall, [numMdls n_dist]);
    
    
    f = figure;
    set(gcf,'Position',get(0,'Screensize'))
    hold on
    for j = 1:numMdls
        plot(1-dec_acc_overall(j, :), 'Color', cols(j, :), 'LineWidth', 1.5)
    end
    plot(1-(1./[2:n_dist+1]), '--k', 'LineWidth', 2)
    lgnd = legend(X);
    % 0.908,0.676,0.089,0.249
    lgnd.Position = [0.908,0.676,0.089,0.249];
    lgnd.FontSize = 16;
    filename = [diskPath filesep taskPath filesep 'DecErr_' num2str(ndim) 'featsLinePlot_' method '_highNC_VC'];
    sgtitle(['Decoding error across models '  num2str(ndim) ' features'], 'FontSize',Fontsize, 'FontWeight', 'bold');
    ylabel('Decoding error', 'FontSize',Fontsize, 'FontWeight', 'bold');
    xlabel('Number of distractors', 'FontSize',Fontsize, 'FontWeight', 'bold')
    set(gca,'FontSize',Fontsize, 'FontWeight', 'bold')
    
    
end

% print(f, filename, '-dpng', '-r0')

% close all



