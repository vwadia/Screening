
setDiskPaths

%% load in all IT cells and sent paths

filePath = 'predictedFeatures';

taskPath = 'Object_Screening';
% load([diskPath filesep 'Object_Screening' filesep 'ITRespCells_500stim_Scrn'])
% load([diskPath filesep 'Object_Screening' filesep 'ITCells_500stim_Scrn_SigRamp'])
% load([diskPath filesep 'Object_Screening' filesep 'ITCells_500stim_Scrn_SigRamp_basicMethod_2Stds'])

load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500stim_Scrn_SigRamp']); merged = 1;

% load([diskPath filesep taskPath filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat']); merged = 0;
% idx = zeros(length(strctCells), 1);
% 
% for cellIndex = 1:length(strctCells)
%     
%     p_info = strctCells(cellIndex).pvalRamp;    
%     
%     idx(cellIndex) = sum(p_info < 0.01);
% end
% 
% strctCells = strctCells(idx == length(p_info));
% responses = responses(idx == length(p_info), :);
% psths = psths(idx == length(p_info), :);

%% load in image params

options.screenType = 'Object';
stimDir = '500Stimuli';
layermat = 'fc6';
% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
params_full = scores{1};
% [coeff, score, ~, ~, ~, mu] = pca(params_full, 'Centered', false);
[coeff, score, ~, ~, ~, mu] = pca(params_full);
params = score(:, 1:50);
python = 0;


% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_' layermat '_' stimDir '.mat']]); matmean = 0;
% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_MatlabMean_' layermat '_' stimDir '.mat']]); matmean = 1;
% [coeff, score, ~, ~, ~, mu] = pca(feat);
% params = score(:, 1:50);
% python = 1;

%% predict features - takes ~70 seconds
tic
% create resp_mat
resp_mat = [];
for idx = 1:length(responses)
    resp_mat(:, idx) = responses{idx, 1};
end

ndim = 50;
all_inds = [1:500]';
% adding column of ones
resp_mat = [resp_mat ones(size(resp_mat, 1), 1)];

for dim = 1:ndim % per dimension
    for idx = 1:length(all_inds) % per image
        
        % using leave one out
        ind_train = setdiff(all_inds, idx);
        Y_train = params(ind_train, dim);
        
        X_train = resp_mat(ind_train, :); % images x cells includes column of 1s
        
        X_test = resp_mat(idx, :);
        
        % solve for coefficients
        [b, ~, ~, ~, stats] = regress(Y_train, X_train);
        
        % predict
        pred_feat(idx, dim) = X_test*b;
        
    end
end
toc
%%
% using pseudo inverse
V_inv = pinv(coeff(:, 1:50));
pred_feat_full = pred_feat*V_inv + repmat(mu, [500 1]);

if python
    if ~matmean
        fname = [diskPath filesep taskPath filesep filePath filesep 'predicted_features_python_pinv_plusmean'];
%         save([diskPath filesep taskPath filesep filePath filesep 'predicted_features_python_pinv_plusmean'], 'pred_feat_full');
    else
        fname = [diskPath filesep taskPath filesep filePath filesep 'predicted_features_python_matmean_pinv_plusmean'];
%         save([diskPath filesep taskPath filesep filePath filesep 'predicted_features_python_matmean_pinv_plusmean'], 'pred_feat_full');
    end
else
    fname = [diskPath filesep taskPath filesep filePath filesep 'predicted_features_matlab_newmean_pinv_plusmean'];
%     save([diskPath filesep taskPath filesep filePath filesep 'predicted_features_matlab_newmean_pinv_plusmean'], 'pred_feat_full');
end

if merged
    fname = [fname '_mergedCells'];
end
save(fname, 'pred_feat_full');


% % using coeffs
% pred_feat_full = pred_feat*coeff(:, 1:50)' + repmat(mu, [500 1]);
% 
% if python
%     save([diskPath filesep taskPath filesep filePath filesep 'predicted_features_python_newmean_multbycoeff_plusmean'], 'pred_feat_full');
% %     save([diskPath filesep taskPath filesep filePath filesep 'predicted_features_python_newmean_multbycoeff'], 'pred_feat_full');
% else
%     save([diskPath filesep taskPath filesep filePath filesep 'predicted_features_matlab_newmean_multbycoeff_plusmean'], 'pred_feat_full');
% %     save([diskPath filesep taskPath filesep filePath filesep 'predicted_features_matlab_newmean_multbycoeff'], 'pred_feat_full');
% end

%% computing normalized distance for GAN

% formula -> norm(fc6_recon - fc6_original)/norm(fc6_bpr - fc6_original)
% fc6_recon = decoded 50D features (pred_feat)
% fc6_original = original features from python (params)
% fc6_bpr = original 4096 features projected into 50D space (have to compute)

load(fname)
% compute fc6_bpr 
proj = params*coeff(:, 1:ndim)' + repmat(mu, [length(all_inds) 1]);


norm_dist = [];
for im = all_inds'
    
    fc6_orig = feat(im, :);
    
    fc6_recon = pred_feat_full(im, :);
    
    fc6_bpr = proj(im, :);
    
    norm_dist(im) = abs(fc6_recon - fc6_orig)/abs(fc6_bpr - fc6_orig);
    
end

%% histogram - see decode_images_largedatabase for a clean way to do this with colors


%% Test decoding accuracy - NN and cosine similarity

metric = 'nn';
% metric = 'cosine_dist';

n_repeats = 1000;
n_comp = 50; 
dec_acc = zeros(size(pred_feat, 1), 1);
% n_distr = 1;

for n_distr = 1:n_comp % takes ~35 min to run when n_comp = 50
    tic
    dec_acc_im = zeros(size(pred_feat, 1), n_repeats);

    for im = 1:size(pred_feat, 1) % for each reconstructed image'
        dist = zeros(1, n_distr+1); 
        for rep = 1:n_repeats
            
            % sample images
            % ensure all images including the target are in sample_ids
            sample_set = setdiff(1:size(params, 1), im);
            sample_ids = [randsample(sample_set, n_distr, false) im];
            
            
            sample_ims = params(sample_ids, :);
            
            % test image
            target_im = pred_feat(im, :);
            
            % compute distances of 50 images to target
            if strcmp(metric, 'nn')
                dist = [];
                dist = vecnorm((sample_ims - target_im)'); % vecnorm finds norm of each column so do it on transpose
            
                % find minimum
                [comp_dist, comp_idx] = min(dist);
                
                if sample_ids(comp_idx) == im
                    dec_acc_im(im, rep) = 1;
                end
            
            elseif strcmp(metric, 'cosine_dist')
                dist = [];
                dist = sample_ims*target_im'; % project target onto samples

                % find maximum
                [comp_dist, comp_idx] = max(dist); % larger proj value = smaller angle between vectors

                if sample_ids(comp_idx) == im
                    dec_acc_im(im, rep) = 1;
                end
            end
            
        end
        dec_acc(im, n_distr) = sum(dec_acc_im(im, :))/n_repeats;
        
    end
    disp(['Finished decoding run with ' num2str(n_distr) ' distractors']);
    toc
end



%% figure - dec acc as a function of distractors

dec_acc_nn = dec_acc;
% dec_acc_cos = dec_acc;
% if ~exist('dec_acc_nn', 'var') || ~exist('dec_acc_cos', 'var')
%     load([diskPath filesep taskPath filesep 'DecAccuracies50Dist_NN_Cos'])
% end
chance = ones(1, n_comp)./[2:n_comp+1];

f = figure;
hold on
plot(mean(dec_acc_nn, 1)', 'LineWidth', 2);
% plot(mean(dec_acc_cos, 1)', 'LineWidth', 2);
plot(chance, '--k', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend({'Nearest Neighbour', 'Cosine Similarity'})
title({'Decoding accuracy all IT cells'});
xlabel('Number of distractors')
ylabel('Decoding accuracy')

% print(f, [diskPath filesep taskPath filesep 'ObjectFeature_Dec_acc_wdistractors_ITCells'], '-dpng', '-r0')


%% figure scatter plot of predicted and actual features (PCs)

% note that this doesn't look great - the model performs badly ona few
% select images makeing the slope of the best fit line is 0.58 before ReLu and then 0.86 after ReLu. 
% potential fix: re-compute axes for all cells with Doris' old school
% latency computation method? --> this doens't make a difference
% EDIT: using the python parameters does though

for pc_num = 1:5

f2 = figure; 
hold on
x = params(:, pc_num);
y = pred_feat(:, pc_num);
if pc_num == 1
    scatter(x, y, 'r', 'filled')
elseif pc_num == 2
    scatter(x, y, 'b', 'filled')
end
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlim([min(x) max(x)])
ylim([min(y) max(y)])
x_lim = double(xlim);
y_lim = double(ylim);
% ylim([-800 500]);
xlabel(['Actual parameter (PC' num2str(pc_num) ')'])
ylabel(['Predicted parameter (PC' num2str(pc_num) ')'])
text(x_lim(1)*0.80, y_lim(2)*0.75, '(n = 500)', 'FontSize', 14,'FontWeight', 'bold');

c = round(corr(x, y), 2);
text(x_lim(1)*0.80, y_lim(2)*0.88, ['Correlation = ' num2str(c)], 'FontSize', 14,'FontWeight', 'bold');

% best fit line
% Get coefficients of a line fit through the data.
coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% plot(xFit, yFit, 'g-', 'LineWidth', 2); % Plot fitted line.
plot(xFit, yFit, 'k-', 'LineWidth', 2); % Plot fitted line.

% filename = [diskPath filesep taskPath filesep 'ModlFeat_ActualFeat_comparison_PC' num2str(pc_num)];
% print(f2, filename, '-dpng', '-r0')
% close all
end

%% Scatter of correlation values of all dimensions - predicted vs normal
f3 = figure; 

c = corr(params, pred_feat); % use this 
% c = corr(params_full, pred_feat_full);
% c = corr(feat, pred_feat_full);

cmap = zeros(size(c, 1), 3);
% cmap(1, :) = [1 0 0]; 
% cmap(2, :) = [0 0 1]; 

scatter(1:size(c, 1), diag(c),'filled', 'CData', cmap); 
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('PC dimension')
ylabel('Correlation value')

% filename = [diskPath filesep taskPath filesep 'Corr_ModelDims_ActualDims'];
% print(f3, filename, '-dpng', '-r0')
% close all

% xl = [-800:450];
% yl = 0.8571*xl;
% 
% plot(xl, yl, 'b-', 'LineWidth', 2) 


%% SCRATCH CODE
%% 

%% figure
% load(fname)
% if merged
%     if matmean % numbers here chosen manually
%         nd_1 = norm_dist(norm_dist < 1.17);
%         nd_3 = norm_dist(norm_dist > 1.26 & norm_dist < 1.41);
%         nd_large = norm_dist(norm_dist > 1.41);
%         
%     else
%         nd_1 = norm_dist(norm_dist < 1.17);
%         nd_3 = norm_dist(norm_dist > 1.26 & norm_dist < 1.41);
%         nd_large = norm_dist(norm_dist > 1.41);
%     end
%     f = figure;
%     hold on
%     histogram(norm_dist , 'FaceColor', [0 0 0])
%     if matmean
%         histogram(nd_1, 3, 'FaceColor', 'r')
%         histogram(nd_3, 'BinEdges', 1.26:0.03:1.41, 'FaceColor', 'b')
%         histogram(nd_large, 'BinEdges', 1.41:0.03:1.89, 'FaceColor', [1 1 1], 'EdgeColor', 'k')
%         
%     else
%         histogram(nd_1, 'BinEdges', 1.05:0.03:1.17, 'FaceColor', 'r')
%         histogram(nd_3, 'BinEdges', 1.26:0.03:1.41, 'FaceColor', 'b')
%         histogram(nd_large, 'BinEdges', 1.41:0.03:1.92, 'FaceColor', [1 1 1], 'EdgeColor', 'k')
%         
%     end
%     xlabel('Normalized distance')
%     ylabel('Number of images')
%     xlim([0.9 2])
%     set(gca, 'FontSize', 14, 'FontWeight', 'bold')
%     if matmean
%         filename = [diskPath filesep taskPath filesep 'Hist_NormDist_GAN_Matmean_SigRampCells'];
%     else
%         filename = [diskPath filesep taskPath filesep 'Hist_NormDist_GAN_SigRampCells'];
%     end
%     if merged
%         filename = [filename '_merged'];
%     end
%     
%     
%     title('Normalized distance of decoded images')
%     
%     
%     print(f, filename, '-dpng', '-r0')
% % end
%% going from 50D --> 4096D (my bootleg method)

% load in full params
% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'allscores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
% params_full = scores{1}; % fc6 before relu
% 
% % params_full = params_full - mean(params_full, 1);
% [coeff, score, ~, ~, explained, mu] = pca(params_full);
% 
% rec = score(:, 1:50)*coeff(:, 1:50)' + repmat(mu, [500 1]);
% recp = params*coeff(:, 1:50)' + repmat(mu, [500 1]);
% 
% V_inv = pinv(coeff(:, 1:50));
% V_inv2 = pinv(coeff);
% [U, S, V_T] = svd(params_full, 'econ');
% 
% V = V_T(:, 1:50)';
%%

% % all give 0.3427 as mean(diag(corr(pred_feat_full, params_full)))
% % pred_feat_full = pred_feat*V_inv + mean(params_full, 1);
% pred_feat_full = pred_feat*V_inv + repmat(mu, [500 1]);
% % pred_feat_full = pred_feat*V_inv2(1:50, :) + mean(params_full, 1);
% % pred_feat_full = pred_feat*coeff(:, 1:50)' + mean(params_full, 1);
% 
% 
% save([diskPath filesep taskPath filesep 'predicted_features_python_pinv_plusmean'], 'pred_feat_full');
