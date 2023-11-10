
% computing STA consistency
% load in all the cells
% Split the responses into 2 halves - use each to compute an axis then correlate them
% repeat this procedure 500x for each cell and take the median
% make barplot (STA consistency vs # of neurons)

% vwadia July 2022


setDiskPaths

taskPath = 'Object_Screening';

layermat = 'fc6';

% stimDir = '500Stimuli';
stimDir = '1593Stimuli';

indiv_sess = 0; % singling out individual session
sess = 'P81CS_AM';

%% P73 session with 1593 images
if strcmp(stimDir, '1593Stimuli')
    
%     load([diskPath filesep taskPath filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328'...
%         filesep 'ITCells_1593Stim_Scrn_basicMethod_2Stds']) 
    load([diskPath filesep taskPath filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328'...
        filesep 'ITCells_1593Stim_Scrn'])
    imageIDs = [1:1593]';
    
    
    load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
    
elseif strcmp(stimDir, '500Stimuli')
    
%     load([diskPath filesep taskPath filesep 'ITCells_500stim_Scrn_SigRamp_basicMethod_2Stds'])
    load([diskPath filesep taskPath filesep 'MergedITCells_500stim_Scrn_SigRamp'])
    imageIDs = [1:500]';
    
    
    
    if indiv_sess
        % single out individual session
        ids = zeros(length(strctCells), 1);
        for cellIndex = l(strctCells)
            
            if strcmp(strctCells(cellIndex).SessionID, sess)
                ids(cellIndex) = 1;
            end
            
            
        end
        
        strctCells(ids == 0) = [];
        responses(ids == 0, :) = [];
        psths(ids == 0, :) = [];
        
%         sess = 'P81CS_1';
        
    else
        sess = 'allCells';
    end
    
    load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
    feat = scores{2};
    
    % load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_MatlabMean_' layermat '_' stimDir '.mat']]);
    [coeff, score, ~, ~, ~, mu] = pca(feat);%, 'Centered', false); % score is much higher range
    params = score(:, 1:50);
end

%%
c_sta = zeros(length(strctCells), 1);
c_lr = zeros(length(strctCells), 1);



% method = 'sta';
% method = 'linear_regression';

n_reps = 100;


ax_1_ws = [];
ax_2_ws = [];
ax_1_lr = [];
ax_2_lr = [];

for cellIndex = 1:length(strctCells)
    
    for rep = 1:n_reps
        hS_1 = sort(randsample(imageIDs, floor(length(imageIDs)/2), false));
        hS_2 = imageIDs(~ismember(imageIDs, hS_1));
        
        resp_1 = responses{cellIndex, 1}(hS_1, 1);
        resp_2 = responses{cellIndex, 1}(hS_2, 1);
        
        parms_1 = params(hS_1, :);
        parms_2 = params(hS_2, :);
        
        
        [ax_1_ws(rep, :), ~] = Utilities.ObjectSpace.analysis_STA(resp_1, parms_1, 'sta');
        [ax_2_ws(rep, :), ~] = Utilities.ObjectSpace.analysis_STA(resp_2, parms_2, 'sta');
        
        [ax_1_lr(rep, :), ~] = Utilities.ObjectSpace.analysis_STA(resp_1, parms_1, 'linear_regression');
        [ax_2_lr(rep, :), ~] = Utilities.ObjectSpace.analysis_STA(resp_2, parms_2, 'linear_regression');
    end
    
    c_sta(cellIndex) = mean(diag(corr(ax_1_ws', ax_2_ws')));
    c_lr(cellIndex) = mean(diag(corr(ax_1_lr', ax_2_lr')));
    disp(['Finished consistency computation for cell no ' num2str(cellIndex)]);
end

ind_con = c_sta > 0.50;
conCells = strctCells(ind_con);
ind_hom = c_sta < 0.30;
homCells = strctCells(ind_hom);

%% figure

f = figure; 
hold on
histogram(c_sta', 10, 'BinEdges', 0.2:0.05:1, 'FaceColor', [0.6350 0.0780 0.1840]);
% histogram(c_lr_res', 10, 'FaceColor', [0.6350 0.0780 0.1840]);
histogram(c_lr', 10, 'BinEdges', 0.2:0.05:1, 'FaceColor', [0.4940 0.1840 0.5560]);
% ylim([0 42])
lgnd = legend({'Weighted Sum', 'Linear Regression'});
xlabel('STA consistency');
ylabel('No of neurons');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% filename = [diskPath filesep taskPath filesep 'STAConsistency_Method_Comparison_MeanValPercell_pythonParams'];
% print(f, filename, '-dpng', '-r300')
% close all

%%
% 
% f = figure; 
% hold on
% histogram(c_sta', 10, 'FaceColor', [0.6350 0.0780 0.1840]);
% % histogram(c_lr_res', 10, 'FaceColor', [0.6350 0.0780 0.1840]);
% histogram(c_lr', 10, 'FaceColor', [0.4940 0.1840 0.5560]);
% ylim([0 42])
% % lgnd = legend({'Weighted Sum', 'Linear Regression'});
% xlabel('STA consistency');
% ylabel('No of neurons');
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');


%% exp vriance fig for AexNet
% [coeff, score, ~, ~, explained, ~] = pca(params_full);
% cs = cumsum(explained);
% ndim = 50;
% f2 = figure; 
% hold on
% scatter(1:250, cs(1:250), 'k')
% plot([ndim ndim], [0 cs(ndim)], '--r', 'LineWidth', 1.5)
% plot([0 ndim], [cs(ndim) cs(ndim)], '--r', 'LineWidth', 1.5)
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% xlabel('Number of PCs')
% ylabel('Explained Variance')
% filename = [diskPath filesep taskPath filesep 'ExpVar_AlexNet_features'];
% print(f2, filename, '-dpng', '-r300')

%% grabbing projection values for cells (rank ordered by FR) - for HSU talk

setDiskPaths

taskPath = 'Object_Screening';

layermat = 'fc6';

stimDir = '500Stimuli';

load([diskPath filesep taskPath filesep 'MergedITCells_500stim_Scrn_SigRamp'])
imageIDs = [1:500]';

load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
options.screenType = 'Object';
options.ind_train = imageIDs;
options.reorder = 1; % IMPORTANT have to set this here
% choos specific cell
% for cellIndex = l(strctCells)
 

% 60/62 = P79 Fingerprint screen - channel 218 cell 1209
% 155/159 = P81 ClosedLoop rescreen - channel 211 cell 658
cellIndex = 63; % 60 = P79 Fingerprint screen - channel 218 cell 1209
projVals = Utilities.ObjectSpace.grab_AxisProj_Values(responses{cellIndex, 1}, params, options);
    
%% Now plot
f3 = figure; 
hold on
% color = [0 0.2 0];

vals_to_plot_top = projVals(end:-1:end-19, 1);
vals_to_plot_bottom = projVals(20:-1:1, 1);

scatter(21:40, flipud(vals_to_plot_top), 30, [1 0 0], 'filled');    
scatter(1:20, flipud(vals_to_plot_bottom), 30, [0 0 1], 'filled') 
% lgnd = legend({'Top stimuli', ' Bottom Stimuli'});
xlabel('Stim ID')
ylabel('Projection Value - Preferred axis')
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
filename = [diskPath filesep taskPath filesep 'ProjVals_vs_FR_Top&BottomStim_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name)];
print(f3, filename, '-dpng', '-r300')


% close all











