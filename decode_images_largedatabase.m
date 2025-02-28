setDiskPaths

taskPath = 'Object_Screening';
% load([diskPath filesep 'Object_Screening' filesep 'ITRespCells_500stim_Scrn'])
load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500stim_Scrn_SigRamp'])

strctCELL = struct2cell(strctCells')';  

%% load in image params
ndim = 50;

options.screenType = 'Object';
stimDir = '500Stimuli';
layermat = 'fc6';

load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
params_full = scores{1}; layer = 'fc6';
% params_full = scores{2}; layer = 'relu6';
% params_full = scores{3}; layer = 'fc7';
% params_full = scores{4}; layer = 'relu7';

% [coeff, score, ~, ~, explained, mu] = pca(params_full);
[coeff, score, ~, ~, ~, mu] = pca(params_full);
params = score(:, 1:ndim);


%%
tic 
% create resp_mat
resp_mat = [];
for idx = 1:length(responses)
    resp_mat(:, idx) = responses{idx, 1};
end

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

%% image decoding - via projection of image set into space create by my ims
%% Does this make sense? Yes but have to remove my images - note that there are very close copies 

%% Finding overlap between both stimulus sets (to remove the 500 images patients saw)

% NOTE there are 2 sources of overlap
% 1) Repeated images in the big set
% 2) the 500 stim images being in the set

% so to remove all images you need to find duplicates AND then check for
% the stim ims
% March2023 vwadia


addpath('ObjectSpace');

% load allimage_resp2.mat % loads im2all - 224x224x19606 image matrix
% load([diskPath filesep taskPath filesep 'LargeImageSet_NoRepeats'])
load([diskPath filesep taskPath filesep 'LargeImageSet_17856ims'])
imageIDs = [1:500]';
% ovrLap_imIDs = [];
% im2all_woStim = im2all;
% im2all_woStim(:, :, ovrLap_imIDs) = [];
% save([diskPath filesep taskPath filesep 'LargeImageSet_NoRepeats'], 'im2all_woStim', '-v7.3')

%% more sophisticated overlap removal 

% a lot of images are visually the same but have a couple of pixels different so don't get removed
% Use SSIM to find the ones super close to stim images (for starters)
bigSetIds = 1:size(im2all_woStim, 3);

ssimVals_imIDs = cell(length(imageIDs), 1);
% pathStim = [diskPath filesep taskPath filesep '500Stimuli'];
% im500 = Utilities.getImageDescriptions(pathStim, 224, 0, 0);
% im500 = uint8(im500); % Leaving as single produces huge artificial differences between images
    
% % finding ims stimilar to stim - took 4 hours with 4 cores
% for ctr = 1:length(imageIDs)
%     targ = im500(:, :, ctr);
%     tic
%     parfor ctr2 = 1:size(im2all_woStim, 3)
%         comp = im2all_woStim(:, :, ctr2);
%         simVal(ctr2) = ssim(targ, comp); 
%     
%     end
%     ssimVals_imIDs{ctr} = simVal;
%     toc
% end

% finding ims similar to other ims in big set
for ctr = 1:length(im2all_woStim)
    targ = im2all_woStim(:, :, ctr);
    tic
    parfor ctr2 = ctr:size(im2all_woStim, 3)
        comp = im2all_woStim(:, :, ctr2);
        simVal(ctr2) = ssim(targ, comp); 
        
    end
    ssimVals_imIDs{ctr} = simVal;
    toc
end


% save([diskPath filesep taskPath filesep 'SSIMValues_500Stimvs17856Set'], 'ssimVals_imIDs', '-v7.3')
%%
closeImIds = {};
% bigSetIds = 1:17856;
parfor ii = 1:length(ssimVals_imIDs)
    nxt = ssimVals_imIDs{ii};
    if sum(nxt > 0.8) > 0
        closeImIds{ii} = bigSetIds(nxt > 0.8);
    end    
end

y = {};
susIds = cellfun(@(x) x', closeImIds, 'UniformOutput', false)';
susIds2 = unique(cell2mat(susIds));

im2all_woStim = uint8(im2all);
im2all_woStim = im2all_woStim(:, :, ~ismember(bigSetIds, susIds2));

%% Removing overlap 
% takes 242s for 18356 images
% reduces set further to 17856 images (March 24th, 2023) 
% im2all = im2all_woStim;
if ~exist('im2all_woStim', 'var')
% im2all = im2all_woStim;
    tic
    ovrLap_imIDs = cell(length(imageIDs), 1);
    ssimVals_imIDs = cell(length(imageIDs), 1);
    
    pathStim = [diskPath filesep taskPath filesep '500Stimuli'];
    im500 = Utilities.getImageDescriptions(pathStim, 224, 0, 0);
    im500 = uint8(im500);
    
    
    for ctr = 1:length(imageIDs)
        
        targ = im500(:, :, ctr);
        rp = 1;
        tic
        for ctr2 = 1:size(im2all, 3)
            
            comp = im2all(:, :, ctr2);
            simVal = ssim(targ, comp);
            
            if isequal(targ, comp) || simVal  > 0.8
                
                ovrLap_imIDs{ctr}(1, rp) = ctr2;
                rp = rp+1;
                
            end
            ssimVals_imIDs{ctr}(1, ctr2) = simVal;
            
        end
        toc
    end
    
    
    repVals = [];
    for val = 1:length(ovrLap_imIDs)
        
        repVals = horzcat(repVals, ovrLap_imIDs{val});
        
    end
    
    ovrLap_imIDs = repVals;
    im2all_woStim = im2all;
    im2all_woStim(:, :, ovrLap_imIDs) = [];
    toc
end


%% compute alexnet features manually for large image set
tic
rem_overlap = 1;

if strcmpi(layer, 'fc6')
    load([diskPath filesep taskPath filesep 'resp_AlexNet_fc6BeforeReLu_17856ims'])
    
    resp_wholeimage = resp_wholeimage( ~ismember(bigSetIds, susIds2), :);
elseif strcmpi(layer, 'relu6')
    if ~rem_overlap
        im2all_woStim = im2all;
    end
    
    grayImages = [];
    for i = 1:size(im2all_woStim, 3)
        
        grayImages(:, :, i) = single(imresize(im2all_woStim(:, :, i), [227 227]));
        
    end
    grayImages = reshape(grayImages, [size(grayImages, 1) size(grayImages, 2) 1 size(grayImages, 3)]);
    grayImages = repmat(single(grayImages), [1 1 3 1]);
    
    
    lexnet = alexnet;
    resp_wholeimage = activations(lexnet, grayImages, layer, 'OutputAs', 'rows');
    
end
proj_into_500 = (resp_wholeimage - repmat(mu, [size(resp_wholeimage, 1) 1]))*coeff; % coeff = PCs of 500 object space
proj = proj_into_500(:, 1:ndim); % now the images are projected into the space built by my 500 ims


toc

% save([diskPath filesep taskPath filesep 'resp_AlexNet_fc6BeforeReLu_largeSetNoRepeats'], 'resp_wholeimage', '-v7.3') 
%% or load in precomputed features

% rem_overlap = 1;
% 
% % load resp_Alexnet.mat % loads responses of Alexnet to 2 conv layers and fc6/7 pre and post relu 
% load([diskPath filesep taskPath filesep 'resp_Alexnet_noRepeats.mat']);
% 
% % response of fc6 for 18356 images
% resp_wholeimage = resp{3};% 3 is for fc6 before relu  
% resp_wholeimage = resp_wholeimage';
% 
% % remove my images
% if rem_overlap
%     resp_wholeimage(ovrLap_imIDs, :) = [];     
% end
% 
% % [coeff2, score2, ~, ~, explained2, mu2] = pca(resp_wholeimage);
% % proj = score2(:, 1:50);
% 
% % projection into space built by stim ims
% proj_into_500 = (resp_wholeimage - repmat(mu, [size(resp_wholeimage, 1) 1]))*coeff; % coeff = PCs of 500 object space
% proj = proj_into_500(:, 1:ndim); % now the images are projected into the space built by my 500 ims

% save([diskPath filesep taskPath filesep 'resp_Alexnet_17856ims.mat'], 'resp', '-v7.3');


%% Find best possible reconstruction - images in the lasrge set with min distance to stim
imageIDs = [1:500];

n_stim = 500;

for im = 1:n_stim
    
    % target image
    target = params(im, :);
    
    % compute distances to all projected images
    dist = proj - repmat(target, [size(proj, 1), 1]);
    
    dd(im, :) = vecnorm(dist');
    % find min
    [d1(im) idx1(im)] = min(dd(im, :));
    
end

if rem_overlap
    ims_bpr = im2all_woStim(:, :, idx1); % select the images;
else
    ims_bpr = im2all(:, :, idx1); % select the images;
end

para_bpr = proj(idx1,:); % parameters of the closest match to all stim images;

%% write out the images for comparison
% 
if ~rem_overlap
    writePath = [diskPath filesep taskPath filesep 'BPRecons500'];
else
    writePath = [diskPath filesep taskPath filesep 'BPRecons500_woStim'];    
end

if ~exist(writePath)
    mkdir(writePath)
end

for i = 1:n_stim
    
    im = ims_bpr(:, :, i);
    
    filename = sprintf('%04d.tif', idx1(i)); % keep the original names
    outputFileName = [writePath filesep [num2str(i) '_' filename]];
    imwrite(im, outputFileName);

end


%% Find decoded images
ddec = [];

for im = 1:n_stim
    
    % target image
    target = pred_feat(im, :);
    
    % compute distances to all projected images
    dist = proj - repmat(target, [size(proj, 1), 1]);
    
    ddec(im, :) = vecnorm(dist');
    
    % find min
    [d2(im) idx2(im)] = min(ddec(im, :));
    
end

if rem_overlap
    ims_decoded = im2all_woStim(:, :, idx2); % select the images;
else
    ims_decoded = im2all(:, :, idx2); % select the images;
end

para_recon = proj(idx2,:); % parameters of the closest match to all stim images;

%% write out the images for comparison

if ~rem_overlap
    writePath = [diskPath filesep taskPath filesep 'DecodedIms500'];
else
    writePath = [diskPath filesep taskPath filesep 'DecodedIms500_woStim'];    
end
if ~exist(writePath)
    mkdir(writePath)
end

for i = 1:n_stim
    
    im = ims_decoded(:, :, i);
    
    filename = sprintf('%04d.tif', idx2(i)); % keep the original names
    outputFileName = [writePath filesep [num2str(i) '_' filename]];
    imwrite(im, outputFileName);

end

%% Computing normalized distance
 
norm_dist = []; %norm(v_recon - v_original)/norm(v_bestpossrecon - v_original)

for im = 1:n_stim
    
    v_orig = params(im, :); % parameters of original image
    
    v_recon = para_recon(im, :);
    
    v_bestPossibleRecon = para_bpr(im, :); % is this right?  
    
    norm_dist(im) = norm(v_recon - v_orig)/norm(v_bestPossibleRecon - v_orig);
    
end

%% plot histogram

sctns = round(prctile(norm_dist, [33, 67]), 1); % 1x2 vector

f = figure;
hold on
histogram(norm_dist(norm_dist <= sctns(1)), 1:0.1:sctns(1), 'EdgeColor', [1 1 1], 'FaceColor', [1 0 0], 'FaceAlpha', 1);
histogram(norm_dist(norm_dist > sctns(1) & norm_dist <= sctns(2)), sctns(1):0.1:sctns(2),'EdgeColor', [1 1 1], 'FaceColor', [0 0 0], 'FaceAlpha', 1);
histogram(norm_dist(norm_dist > sctns(2)), sctns(2):0.1:7,'EdgeColor', [1 1 1], 'FaceColor', [0 0 1], 'FaceAlpha', 1);
% histogram(norm_dist(norm_dist > sctns(2)), sctns(2):0.1:4,'EdgeColor', [1 1 1], 'FaceColor', [0 0 1], 'FaceAlpha', 1);

% without white edges so you can show full dist
% histogram(norm_dist(norm_dist <= sctns(1)), 1:0.1:sctns(1), 'FaceColor', [1 0 0], 'FaceAlpha', 0.8);
% histogram(norm_dist(norm_dist > sctns(1) & norm_dist <= sctns(2)), sctns(1):0.1:sctns(2), 'FaceColor', [0 0 0], 'FaceAlpha', 0.8);
% histogram(norm_dist(norm_dist > sctns(2)), sctns(2):0.1:15, 'FaceColor', [0 0 1], 'FaceAlpha', 0.8);


% histogram(norm_dist);
% histogram(norm_dist_largedat, 0:0.1:3)
title('Normalized distance of decoded images')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Normalized distance')
ylabel('Number of Images')
x_lim = xlim;
xlim([1 x_lim(2)])

% print(f, [diskPath filesep taskPath filesep 'Hist_NormDist_LargeImageSet_SigRampcells_merged_' layer], '-dsvg', '-r300')


% %% patchShow with decoded images (make this better)
% 
% pathStim = [diskPath filesep taskPath filesep '500Stimuli'];
% stimDir = dir(fullfile(pathStim));
% stimDir = stimDir(~ismember({stimDir.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));
% 
% % gets rid of weird copies
% stimNames = struct2cell(stimDir);
% stimNames = stimNames(1, :)';
% goodStim = ~startsWith(stimNames, '._', 'IgnoreCase', true);
% stimDir = stimDir(goodStim);
% 
% [~, natIdx] = natsortfiles({stimDir.name});
% 
% for i = 1:length(stimDir)
%     
% %     stim(:, :, i) = imread([stimDir(natIdx(image)).folder filesep stimDir(natIdx(image)).name]);
%     stim(:, :, i) = imread([stimDir(i).folder filesep stimDir(i).name]);
%     
% end
% patchFig = Utilities.Plotting.patchShow(stim, ims, 20);
% if rem_overlap
%     filename = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'patched_largeImageSetwoStim_decoding'];
% else
%     filename = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'patched_largeImageSet_decoding'];
% end
% print(patchFig, filename, '-dsvg', '-r300')
% close all



%% SCRATCH CODE
 %% writing out images to look at them
% im2all = im2all_woStim;
BigSetNum = size(im2all_woStim, 3);
writePath = [diskPath filesep taskPath filesep num2str(BigSetNum) 'Stimuli'];
if ~exist(writePath)
    mkdir(writePath)
end

% n_stim = size(im2all, 3);
% idxes = setdiff(1:size(im2all, 3), ovrLap_imIDs);

for i = 1:size(im2all_woStim, 3)
    
    im = im2all_woStim(:, :, i);
    
    filename = sprintf('%04d.tif', i);
    outputFileName = [writePath filesep filename];
    imwrite(im, outputFileName);
    
end


%% remove duplicates 
% im2all = im2all_woStim;
% tic
% dup_ids = zeros(size(im2all, 3), 1);
% 
% for i = 1:size(im2all, 3)
%     
%     targ = im2all(:, :, i);
%     for j = (i+1):length(dup_ids)
%         
%         if isequal(im2all(:, :, j), targ)
%             dup_ids(j) = 1;
%             disp([i, j])
%         end
%         
%     end
%     
% end
% 
% % remove duplicate indices
% im2all(:, :, find(dup_ids == 1)) = [];
% % save([diskPath filesep taskPath filesep 'LargeImageSet_NoRepeats'], 'im2all', '-v7.3')
% % save([diskPath filesep taskPath filesep 'LargeImageSet_RepeatedIndices'], 'dup_ids');
% 
% toc




