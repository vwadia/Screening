%% Script to decode image features based on population responses
% 
%% For a given cell:
% Resp = a*X + b
% where:
% Resp = n_img x 1 vector of responses
% X = para*sta'
%     para = n_img x n_dim matrix of image features
%     sta = mean_sub_resp'*para (preferred axis)

% In this situation we compute the STA using n-1 images, project the full params matrix onto it
% take the training indices of that projection, and the training indices of the response and solve linear equation
% plug coefficients back in to compute test response of test image.

%% To decode image features we want to use population activity

% feat = pop_resp*c + d
% 
% pop_resp = n_img x n_cells responses of all cells to images
% feat     = para = n_img x ndim, features of all images
%    c     = n_cells x ndim matrix of beta  regression coefficients  
%    d     = intercept
% use LOO to fit c, d and predict feature vec of 1 image

% Per feature value (per dimension)
%     use n-1_img x 1 feature vec (y_train), n-1_img x n_cells resp matrix (x_train) 
%     to fit c and d 
%     then predict 


if strcmp(taskStruct.subID, 'P73CS_ParamObj') 
%     tic
%     ndim = 50;
%     load allimage_resp2.mat % produces im2all
%     load resp_Alexnet.mat
%     % response of fc6 for 19606 images
%     layer_resp = resp{3}';% 3 is for fc6 before relu
%     load r_1593
%     [coeff1593,score,latent,tsquared,explained,mu1593] = pca(r_1593);
%     params = score(:, 1:ndim);
%     toc
    load parameters_1593_objects.mat % creates params
elseif strcmp(taskStruct.subID, 'P73CS_ParamObj_500')  
    options.screenType = 'Object';
    stimDir = '500Stimuli';
    layermat = 'fc6';
    load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
elseif strcmp(taskStruct.subID, 'P73CS_Full')
    % load params
    load parameters_2k_synthetic_faces.mat
end

%%

all_inds = screeningData.imageIDs;
ndim = 50;
resp_mat = [];
ctr = 1;
for i = l(screeningData.responses)
    if ~isempty(screeningData.responses{i, 1})
        resp_mat(:, ctr) = screeningData.responses{i, 1};
        ctr = ctr + 1;
    end 
end

%%
% adding column of ones
resp_mat = [resp_mat ones(size(resp_mat, 1), 1)];

for dim = 1:ndim % per dimension
    for idx = l(all_inds) % per image
        
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

%% Find the target images

% destPath =[basePath filesep 'rasters' filesep 'Reconstructions'];
% if ~exist(destPath)
%     mkdir([destPath])
% end
% 
% 
% if strcmp(taskStruct.subID, 'P73CS_ParamObj')
% 
%     generate_objects_50d_obj_space(destPath, layer_resp, coeff1593, mu1593, im2all, ndim, 'recon', pred_feat); 
%     
% elseif strcmp(taskStruct.subID, 'P73CS_Full')
%         std_params = std(params);
% %     for im = 1:size(pred_feat, 1)
%         generate_face_50d_face_space(destPath,  'recon', pred_feat, std_params)
% %     end
% end

%% decoding accuracy

% draw 50 faces randomly out of the params
% compute euclidean distance between all vectors and the target
% rank order them and see if top index is the same as target (yes 1, no 0)
% repeat this procedure 1000 times and compute accuracy

n_repeats = 1000;
n_comp = 50;
dec_acc = zeros(size(pred_feat, 1), 1);
dec_acc_im = zeros(size(pred_feat, 1), n_repeats);
n_distr = 50;
% for n_distr = 1:n_comp-1
    
    for im = 1:size(pred_feat, 1) % for each reconstructed image'
        
        for rep = 1:n_repeats
            
            % sample images
            % ensure all images including the target are in there once
            sample_set = setdiff(1:size(params, 1), im);
            sample_ids = [randsample(sample_set, n_distr, false) im];
            
            
            sample_ims = params(sample_ids, :);
            
            % test image
            target_im = pred_feat(im, :);
            
            % compute distances of 50 images to target
            dist = [];
            for run = 1:n_distr
                dist(run) = norm(sample_ims(run, :) - target_im);
            end
            % find minimum
            [min_dist, min_idx] = min(dist);
            
            if sample_ids(min_idx) == im
                dec_acc_im(im, rep) = 1;
            end
        end
        dec_acc(im) = sum(dec_acc_im(im, :))/n_repeats;
        
    end

% end












