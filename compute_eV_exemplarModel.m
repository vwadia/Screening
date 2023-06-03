% compute_ev_exemplar model


% making gradient descent happen for exemplar model

setDiskPaths

% screenType = 'Face';
screenType = 'Object';

if strcmp(screenType, 'Object')
    taskPath = 'Object_Screening';
    addpath(genpath('ObjectSpace'))
    
    load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat']); RespCells = 1; layermat = 'fc6'; layerpy = 'IT_output'; stimDir = '500Stimuli';
    %     load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500Stim_Scrn_SigRamp.mat']); SigCells = 1; layermat = 'fc6'; layerpy = 'IT_output'; stimDir = '500Stimuli';
    
    outPath = [diskPath filesep taskPath filesep 'ExemplarModel_EV'];
    if ~exist(outPath)
        mkdir(outPath)
    end
    
    imageIDs = [1:500]';
    
    
    load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
    
    % set up order (for noise ceiling and exp var computation)
    stimDir = '500Stimuli';
    imageIDs = [1:500]';
    order1 = repelem(imageIDs, 6);
    order2 = repelem(imageIDs, 4);
    catOrd = {};
    
    % afternoon session that got fucked
    a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];
    a_strct = load([a_basePath filesep 'PsthandResponses']);
    order3 = a_strct.order;
    clearvars a_strct
    
    anovaType = 'CategoryObject';
    faceInds = 134:210;
    objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
    textInds = [264:282 290 400:408];
    vegInds = [211:235 357:399];
    animInds = [1:84 256:263];
    
    for i = 1:3
        
        if i == 1
            order = order1; % 3000 elements
        elseif i == 2
            order = order2; % 2000 elements
        elseif i == 3
            order = order3;
        end
        
        catOrder = zeros(length(order), 1);
        catOrder(ismember(order, faceInds)) = 1;
        catOrder(ismember(order, textInds)) = 2;
        catOrder(ismember(order, vegInds)) = 3;
        catOrder(ismember(order, animInds)) = 4;
        catOrder(ismember(order, objInds)) = 5;
        
        catOrd{i} = catOrder;
    end
    
    
elseif strcmp(screenType, 'Face')
    % second attempt with just SA faces - does this model work at all?
    
    taskPath = 'Face_Screening';
    outPath = [diskPath filesep taskPath filesep 'ExemplarModel_EV'];
    if ~exist(outPath)
        mkdir(outPath)
    end
    
    imageIDs = [1:667]';
    
    
    load([diskPath filesep 'ObjectSpace' filesep 'parameters_2k_synthetic_faces.mat']);
    params = params(1:667, :);
    
    sessPath = ['P73CS' filesep '2000FaceScreen_Session_1_20210326'];
    
    basePath = [diskPath filesep taskPath filesep sessPath];
    
    load([basePath filesep 'PsthandResponses'])
    load([basePath filesep 'strctCells'])
    
    %     load([basePath filesep 'realOrder_P73CS_FullSynthFaces.mat']);
    %     order = imageTexture - 10; % image textures are labelled by the picture num + 10
    %     order = order';
end

%%

tune_e_first = false; % imp for SA faces but really screws up ev in obj space
norm_para = false;

set_norm = true; % imp for objects too - otherwise ev has ~1e10 magnitude
max_norm = 3;

%%
% choose subset of cells to test on
% 0.1672 0.1275 0.2144 0.0579 0.1132 0.0975 - respective ev w/o norm
% respcellSub = [14 101 178 232 277 349]'; % only ones with poritive ev out of first 52

% respcellSub = randsample(length(strctCells), 40); % only ones with poritive ev out of first 52
% cellSubIdx = ismember(1:length(strctCells), respcellSub);


% respcellSub = randsample(length(strctCells), 5);
% cellSubIdx = ismember(1:length(strctCells), respcellSub);
% responses(~cellSubIdx, :) = [];
% psths(~cellSubIdx, :) = [];
% strctCells(~cellSubIdx) = [];
%

%% Modified implementation of Steven's code
%% Without normalizing for all responsive cells this took 14+ hours on 8 cores...

if norm_para
    amp_dim = sqrt(sum(params.^2)); % finding the norm of each dimension 1xndim
    amp_dim(amp_dim == 0) = 1;
    if strcmp(screenType, 'Object')
        params = param_normalize_per_dim(params, amp_dim, length(imageIDs));
    elseif strcmp(screenType, 'Face')
        params = param_normalize(params, amp_dim, length(imageIDs));
    end
end

opts = optimoptions('lsqcurvefit');
% opts.MaxFunctionEvaluations = 2.5e4;
opts = optimset('Display', 'off');

tic
all_inds = imageIDs';
ndim  = size(params, 2);
mn_ctr = 1;
for max_norm = 10
    for cellIndex = 1:length(strctCells)
        
        if ~isnan(responses{cellIndex, 1})
            % grab responses
            obs = responses{cellIndex, 1};
            pred = [];
            parfor idx = all_inds
                %     for idx = all_inds
                
                % train and test indices
                ind_train = setdiff(all_inds, idx);
                ind_test = idx;
                
                % set up training and testing features
                para_train = double(params(ind_train, :));
                para_test = double(params(idx, :));
                
                obs_train = obs(ind_train);
                obs_test = obs(idx);
                
                % compute sta with held out set
                sta = para_train' * obs_train;
                if set_norm
                    sta = sta/norm(sta)*max_norm;
                end
                
                
                if tune_e_first
                    % compute distances between sta and training params
                    dist = vecnorm((para_train - sta'), 2, 2);
                    
                    % regress onto responses for initial guess
                    [beta] = regress(obs_train, [dist ones(length(dist), 1)]);
                    
                    % linear distance needs to take in only 2 parameters
                    exem_beta = [sta; beta];
                    
                    % Initialize exemplar bounds
                    % using ndim so
                    ub = Inf(1, length(exem_beta));
                    %             lb = -Inf(1, length(exem_beta));
                    lb = zeros(1, length(exem_beta));
                    
                    if set_norm
                        exem_beta(end+1) = max_norm;
                        lb(end+1) = 0;
                        ub(end+1) = max_norm;
                    end
                    
                    exem = lsqcurvefit(@linear_distance, exem_beta', para_train, obs_train', lb, ub, opts);
                else
                    exem = sta';
                end
                
                % recompute distances
                dist = vecnorm((para_train - exem(1:ndim)), 2, 2);
                
                % fit 3rd order poly
                c = polyfit(dist, obs_train, 3);
                
                % evaluate polynomial with least squares minimization
                % reset the bounds
                poly_exem = [exem(1:ndim) c];
                ub = Inf(1, length(poly_exem));
                %         lb = -Inf(1, length(poly_exem));
                lb = zeros(1, length(poly_exem));
                
                if set_norm
                    poly_exem(end+1) = max_norm;
                    lb(end+1) = 0;
                    ub(end+1) = max_norm;
                    
                end
                exem_final = lsqcurvefit(@poly_distance, poly_exem, para_train, obs_train', lb, ub, opts);
                
                % use these to compute ev
                d_test = norm(exem_final(1:ndim) - para_test);
                pred(idx, 1) = polyval(c, d_test);
                
                % should save exem_final here?
                
            end
            all_pred(:, cellIndex) = pred;
            error_regress = sum((pred - obs).^2);
            error_total = sum((obs - mean(obs)).^2);
            ev(mn_ctr, cellIndex) = 1 - error_regress./error_total;
            %     ev(1, cellIndex) = 1 - error_regress./error_total;
            
            % define order
            if strcmp(screenType, 'Object')
                
                if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
                    labels = catOrd{1};
                    timelimits = [-0.17, 0.33];
                    stimOffDur = 166.6250;
                    stimDur = ceil(166.6250);
                    order = order1;
                elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
                    labels = catOrd{3};
                    timelimits = [-0.17, 0.33];
                    stimOffDur = 133.4680;
                    stimDur = ceil(266.6250);
                    order = order3;
                else
                    labels = catOrd{2};
                    timelimits = [-0.17, 0.53]; %changed for all cells Oct2022 vwadia
                    stimOffDur = 133.4680;
                    stimDur = ceil(266.6250);
                    order = order2;
                end
                expV(mn_ctr, cellIndex) = Utilities.computeExplainableVariance(psths(cellIndex, :), order, floor(responses{cellIndex, 2}), stimDur, -timelimits(1)*1e3);
                
            else
                expV(mn_ctr, cellIndex) = Utilities.computeExplainableVariance(psths(cellIndex, :), order, floor(responses{cellIndex, 2}), 267, 170);
            end
            
        end
        
        
    end
    
    adj_ev(:, mn_ctr) = ev(mn_ctr, :)./expV(mn_ctr, :);
    
    
    mn_ctr = mn_ctr+1;
    
end
toc

% save([outPath filesep 'RampCells_Norm10_lbInf_ev'], 'ev');


%% helpers
% normalization
function param = param_normalize_per_dim(param, amp_dim, NIMAGE)
%% normalize each dimension separately
param = param./repmat(amp_dim, [NIMAGE 1]);
end

function param = param_normalize(param, amp_dim, ndim1)
%% normalize shape/appearance separately while keeping the relative amplitude within shape or appearance dimensions
%% stevens way - in the cell paper
ndim = size(param, 2);
% para = para./repmat(amp_dim, [NIMAGE 1]);

param(:,1:ndim1)=param(:,1:ndim1) / sqrt(sum(amp_dim(1:ndim1).^2)) / sqrt(2);
param(:, ndim1+1:ndim)=param(:, ndim1+1:ndim) / sqrt(sum(amp_dim(ndim1+1:ndim).^2)) / sqrt(2);
end

% linear distance
function F = linear_distance(exem_betas, features)
% take in exemplar and training features
% compute linear function of distance between features and exem
% INPUTS:
%     1. exemplar and betas as 1 vector - ndim + 2 vector
%     2. features - training ims x ndim params
%
% OUTPUTS:
%     1. updated exemplar which is a function of the linear distance
%     bw exemplar and features
% vwadia May2023

if length(exem_betas) == size(features, 2)+3
    betas(1) = exem_betas(end-2);
    betas(2) = exem_betas(end-1);
    exem = exem_betas(1:end-3);
    des_norm = exem_betas(end);
    
    % normalize if available
    exem = exem./norm(exem)*des_norm;
else
    betas(1) = exem_betas(end-1);
    betas(2) = exem_betas(end);
    exem = exem_betas(1:end-2);
end

for i = 1:size(features, 1)
    F(i) = betas(1) + betas(2)*norm(exem - features(i, :));
end
end

% polynomial distance
function F = poly_distance(exem, features)
% take in exemplar, training features, and poly coefficients
% evaluate function of distance between features and exem
% INPUTS:
%     1. exemplar - ndim vector
%     2. features - training ims x ndim params
%     3. coeffs - coefficients of 3rd degree polyfit
%
% OUTPUTS:
%     1. updated exemplar which is a function of the linear distance
%     bw exemplar and features
% vwadia May2023
if length(exem) == size(features, 2)+5
    c = exem(end-4:end-1);
    des_norm = exem(end);
    d = vecnorm(exem(1:end-5)' - features');
else
    c = exem(end-3:end);
    d = vecnorm(exem(1:end-4)' - features');
end
F = c(1)*d.^3 + c(2)*d.^2 + c(1)*d + c(4);

end




%% My old implementation -  doesn't really work
% n_iter = 1000; % for grad descent
%
% ndim = size(params, 2);
% alpha_e = 0.001; % learning rate for e
% alpha_c0 = 1e-10; % learning rate for coeffs
% alpha_c1 = 1e-8; % learning rate for coeffs
% alpha_c2 = 1e-6; % learning rate for coeffs
% alpha_c3 = 1e-4; % learning rate for coeffs
% tol = 1e-5;
%
% cellExem = nan(length(strctCells), ndim);
% cellCoeffs = nan(length(strctCells), 4);
%
% % per cell
% for cellIndex = 1:length(strctCells)
%
%     % initialize exemplar
%     e = (rand(size(params, 2), 1).*(-800) + 400)'; % random numbers b/w -1 and 1
%     e = e./norm(e);
% %     e = mean(abs(params), 1); % random numbers
%     oldNorm(cellIndex, :) = norm(e);
%
%     % grab responses
%     obs = responses{cellIndex, 1};
%
%     % compute dist
%     dist = vecnorm((params - e), 2, 2);
%
%     % solve for coefficients initially
%     c = polyfit(dist, obs, 3);
%     c0 = c(1); c1 = c(2); c2 = c(3); c3 = c(4);
%
%     stopC0 = false;
%     stopC1 = false;
%     stopC2 = false;
%     stopC3 = false;
% %     tic
%     for iter = 1:n_iter
%
%         dist = vecnorm((params - e), 2, 2);
%
%         % compute predicted responses - evaluate polynomial
%         c = polyfit(dist, obs, 3);
%         pred = polyval(c, dist);
%         c0 = c(1); c1 = c(2); c2 = c(3); c3 = c(4);
%
% %         c0 = c(1); c1 = c(2); c2 = c(3); c3 = c(4);
%
%         % cost function - regression error
%         error_regress = sum((pred-obs).^2);
%         cost = error_regress;
%         delta = pred - obs;
%
%         % gradient wrt e, summing contributions over all stimuli
%         % grad_e = d(cost)/d(pred) * d(pred)/d(dist) * d(dist)/d(e)
%         for i = 1:length(e)
% %             grad_e(i) = sum(2*(delta)) .* sum((3*c0*dist.^2 + 2*c1*dist + c2)) .* sum((e(i) - params(:, i))./dist);
%             grad_e(i) = sum(2*(delta) .* (3*c0*dist.^2 + 2*c1*dist + c2) .* (e(i) - params(:, i))./dist);
%         end
%
%
%         % gradients wrt other polynomial coefficients
%         % these blow up because dist is very large...
% %         grad_c0 = sum(2*delta .* dist.^3);
% %         grad_c1 = sum(2*delta .* dist .^2);
% %         grad_c2 = sum(2*delta .* dist);
% %         grad_c3 = sum(2*delta);
%
%         % update exemplar
%         e_new = e - alpha_e * grad_e;
%
%         % update coeffs
% %         if ~stopC0
% %             c0_new = c0 - alpha_c0 * grad_c0;
% %         end
% %         if ~stopC1
% %             c1_new = c1 - alpha_c1 * grad_c1;
% %         end
% %         if ~stopC2
% %             c2_new = c2 - alpha_c2 * grad_c2;
% %         end
% %         if ~stopC3
% %             c3_new = c3 - alpha_c3 * grad_c3;
% %         end
%
%
%         % check for convergence
%         if max(abs(e_new - e)) < tol
%             break
% %         elseif abs(c0_new - c0) < tol
% %             stopC0 = true;
% %         elseif abs(c1_new - c1) < tol
% %             stopC1 = true;
% %         elseif abs(c2_new - c2) < tol
% %             stopC2 = true;
% %         elseif abs(c3_new - c3) < tol
% %             stopC3 = true;
%         end
%
%         % reassign for next iteration
%          e = e_new;
% %          c0 = c0_new;
% %          c1 = c1_new;
% %          c2 = c2_new;
% %          c3 = c3_new;
%
%     end
% %     disp(iter)
% %     toc
%     cellExem(cellIndex, :) = e;
%     cellCoeffs(cellIndex, :) = [c0 c1 c2 c3];
%
%     % now compute explained variance
%     dist = vecnorm((params - e), 2, 2);
%     c = [c0 c1 c2 c3];
%     pred = polyval(c, dist);
%     obs = responses{cellIndex, 1};
%
%     error_regress = sum( (pred - obs).^2 );
%     error_total = sum((obs-mean(obs)).^2);
%     ev(cellIndex) = 1-error_regress./error_total;
%
% end
