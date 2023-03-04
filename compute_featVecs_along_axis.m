

% Script to compute feature vectors n steps before and after the STA along the axis
% This is to produce stimuli along the axis to see what features cells are responding to
% (Emily's SURF project)
%
% This needs checking
% vwadia July2022

%% load in all IT cells and sent paths

% setDiskPaths
% 
% taskPath = 'Object_Screening';
% load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500stim_Scrn_SigRamp']); splitByCat = 0;
% % load([diskPath filesep 'Object_Screening' filesep 'ITResponses_SigRamp_SplitByCategory_500Stimuli']); splitByCat = 1;
% 
% compress_expand = 0; % compress via PCA and expand
% 
% indiv_sess = 1; % singling out individual session
% sess = 'P81CS';
% 
% if indiv_sess
%     % single out individual session
%     ids = zeros(length(strctCells), 1);
%     for cellIndex = 1:length(strctCells)
% 
%         if strcmp(strctCells(cellIndex).SessionID, sess)
%             ids(cellIndex) = 1;
%         end
% 
% 
%     end
% 
%     strctCells(ids == 0) = [];
%     responses(ids == 0, :) = [];
%     psths(ids == 0, :) = [];
% 
% %     sess = 'P81CS_AM';
% 
% else
%     sess = 'allCells';
% end

%% load in a single session (if it hasn't been merged with full cell set yet)

setDiskPaths
indiv_sess = 1;
compress_expand = 1;
splitByCat = 0;
    
taskPath = 'Object_Screening';
% sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'FingerprintScreening_Session_1_20230111'];
sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115'];

load([sessPath filesep 'PsthandResponses']);
load([sessPath filesep 'strctCells']);
sess = strctCells(1).SessionID;
%% load in params - adjust this for screen type

options.screenType = 'Object';
stimDir = '500Stimuli';
layermat = 'fc6';

% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
% params_full = scores{2}; %

% what are differences between these 2?
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_MatlabMean_' layermat '_' stimDir '.mat']]);
% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_' layermat '_' stimDir '.mat']]);

if compress_expand
    [coeff, score, ~, ~, ~, mu] = pca(feat);%, 'Centered', false); % score is much higher range
    params = score(:, 1:50);
    V_inv = pinv(coeff(:, 1:50));
else
    params = feat;
end
% params = feat;
python = 1;


%% compute feature vectors

% compute STA - there is a 1e5 scale difference between these two.
% sta is bigger, linear_regression is smaller
% should not be an issue if one normalizes - makes sense ans regression
% returns a unit vector basically so if you normalize sta then voila
all_norms = [];

type = '';
% params = params_full; type = '';
% params = params - mean(params); type = '_meanSub_';
% params = recon_params; % is this done correctly??

method = 'sta';
% method = 'linear_regression';

addNoise = 0;

featMat = {};
n_steps = 5;
scale = 1;
n_steps_ortho = 5;

% image_spread = 'PrefOnly'; 
image_spread = 'PrefandOrtho'; 
% image_spread = 'Grid'; 

if strcmp(image_spread, 'PrefOnly')
    stepRangeOrtho = setdiff([-n_steps_ortho:1:n_steps_ortho], 0)*scale;
    stepRange = setdiff([-n_steps:1:n_steps], 0)*scale;
    numcols = 1;
elseif strcmp(image_spread, 'PrefandOrtho')
    stepRangeOrtho = setdiff([-n_steps_ortho:1:n_steps_ortho], 0)*scale;
    stepRange = setdiff([-n_steps:1:n_steps], 0)*scale;
    numcols = 2;
elseif strcmp(image_spread, 'Grid')
    % need the 0 here because preferred axis = 0 ortho coord
    stepRangeOrtho = [-n_steps_ortho:1:n_steps_ortho]*scale;
    stepRange = setdiff([-n_steps:1:n_steps], 0)*scale;
    numcols = length(stepRangeOrtho)+1; % add one column to include ortho axis itself
end


if splitByCat
    catRange = 1:5;
else
    catRange = 0;
end

for cat = catRange
    
    switch cat
        case 0
            group = '';
        case 1
            group = 'faceCells';
            strctCells = faceCells;
            responses = faceCells;
        case 2
            group = 'textCells';
            strctCells = textCells;
            responses = textCells;
        case 3
            group = 'vegCells';
            strctCells = vegCells;
            responses = vegCells;
        case 4
            group = 'animCells';
            strctCells = animCells;
            responses = animCells;
        case 5
            group = 'objCells';
            strctCells = objCells;
            responses = objCells;
    end
    
%     cellIndex = 1;
    for cellIndex = 3 %1:length(strctCells)
        
        if ~isempty(responses{cellIndex, 2})
            resp = responses{cellIndex, 1};
            [sta, ~] = Utilities.ObjectSpace.analysis_STA(resp, params, method);
            
            if ~strcmp(image_spread, 'PrefOnly')
                % Normalize params - note that this is done inside
                % analysis_STA for pref axis
                amp_dim = sqrt(sum(params.^2));
                amp_dim(amp_dim == 0) = 1;
                para = param_normalize_per_dim(params, amp_dim, length(resp));
                
                % note this is different - using non-normalized params
                value_sta_prj = (sta/norm(sta))*params';
                
                para_sub_sta = zeros(size(para));
                for k=1:size(para,1)
                    param_sta_prj = sta*(para(k,:)*sta')/(sta*sta'); % vector of params pojected onto STA
                    para_sub_sta(k,:) = para(k,:) - param_sta_prj; % subtract STA component from param
                end

                % order orthogonal axes by variance explained
                % note norm of COEFF(:,1) = 1
               COEFF = pca(para_sub_sta);
                
                % choose longest
                orthAx = COEFF(:,1)';
            end
            
%             if compress_expand
%                 std_params = std(params);
%             else
                std_params = std(params);
%             end
            
            if strcmp(method, 'linear_regression')
                if strcmp(options.screenType, 'Face') % shape appearance screen
                    
                    % shape-appearance dimensions are flipped between mine and Stevens code
                    % why is this only for lin-reg and not general sta?
                    sta = fliplr(sta)./norm(sta); % normalize length
                    
                else
%                     if compress_expand
%                         sta = sta*V_inv + mu;
%                     end
                    sta = sta./norm(sta);
                    
                end
                all_sta_linReg(cellIndex, :) = sta;
                
            elseif strcmp(method, 'sta')
%                 if compress_expand
%                     sta = sta*V_inv + mu;
%                 end
                sta = sta./norm(sta); % normalize length
                all_sta(cellIndex, :) = sta;
                
            end
            
            if exist('orthAx', 'var')
                all_orthAx(cellIndex, :) = orthAx;
            end
            
            
            if addNoise
                % add noise - optional
                sigmas = 0.05 * sta; % 0.05 = noise level 5%.
                randomNoise = randn(1, length(sta)) .* sigmas;
                % Add noise to a to make an output column vector.
                sta = sta + randomNoise;
            end
            
            
            % create cell array - row of cells per neuron
            % each cell - matrix with steps along (pref axis, nth step along orthogonal axis)
            % Eg. for PrefOnly this matrix will be cells x 2, for
            % PrefandOrtho it will be cells x 3
            % for Grid it will be cells x number of steps along ortho+1
            % the extra cell in each case is for the neuronID
            realStep = 1;
            if strcmp(image_spread, 'PrefOnly')
                for step = stepRange
%                     if step == 0
%                         featMat{cellIndex, 1}(realStep, :) = -sta;
%                         featMat{cellIndex, 1}(realStep+1, :) = sta;
%                         realStep = realStep + 2;
%                     else
                        % if step = 1 will it actually just return STA? Yes that is the actual STA stimulus
                        % the variable 'sta' is a unit vector in that direction
                        featVec = sta.*(std_params*step);
                        if compress_expand
                            featMat{cellIndex, 1}(realStep, :) = featVec*V_inv + mu;
                        else
                            featMat{cellIndex, 1}(realStep, :) = featVec;
                        end
                        realStep = realStep + 1;
%                     end
                end
            elseif strcmp(image_spread, 'PrefandOrtho')
                
                % preferred axis
                realStep = 1;
                for step = stepRange
%                     if step == 0
%                         featMat{cellIndex, 1}(realStep, :) = -sta;
%                         featMat{cellIndex, 1}(realStep+1, :) = sta;
%                         realStep = realStep + 2;
%                     else
                        % if step = 1 will it actually just return STA?
                        % unclear if sta is unit vec * std dev or the raw
                        % sta returned by function
                        featVec = sta.*(std_params*step);
                        if compress_expand
                            featMat{cellIndex, 1}(realStep, :) = featVec*V_inv + mu;
                        else
                            featMat{cellIndex, 1}(realStep, :) = featVec;
                        end
                        
                        
                        realStep = realStep + 1;
%                     end
                end
                
                % othogonal axis
                realStepOrtho = 1;
                for orthostep = stepRangeOrtho
%                     if orthostep == 0
%                         featMat{cellIndex, 2}(realStepOrtho, :) = -orthAx;
%                         featMat{cellIndex, 2}(realStepOrtho+1, :) = orthAx;
%                         realStepOrtho = realStepOrtho + 2;
%                     else
                        featVec = orthAx.*(std_params*orthostep);
                        if compress_expand
                            featMat{cellIndex, 2}(realStepOrtho, :) = featVec*V_inv + mu;
                        else
                            featMat{cellIndex, 2}(realStepOrtho, :) = orthAx.*(std_params*orthostep);
                        end
                                             
                        realStepOrtho = realStepOrtho + 1;
%                     end
                end
                
                
            elseif strcmp(image_spread, 'Grid')
                
                % step along othogonal axis
                realStepOrtho = 1;
                for orthostep = stepRangeOrtho
                    toAdd = orthAx.*(std_params*orthostep);
                    
                    % make a series of steps along pref for each otho step
                    realStep = 1;
                    for step = stepRange
%                         if step == 0 % these are just unit vecs no need
%                         to add
%                             featMat{cellIndex, realStepOrtho}(realStep, :) = -sta + toAdd;
%                             featMat{cellIndex, realStepOrtho}(realStep+1, :) = sta + toAdd;
%                             realStep = realStep + 2;
%                         else                          
                            featVec = sta.*(std_params*step) + toAdd;
                            if compress_expand
                                featMat{cellIndex, realStepOrtho}(realStep, :) = featVec*V_inv + mu;
                            else
                                featMat{cellIndex, realStepOrtho}(realStep, :) = featVec;
                            end
                            
                            realStep = realStep + 1;
%                         end
                    end
                    
                    % move to the next point along ortho ax
                    realStepOrtho = realStepOrtho + 1;
                    
                end
                
                % grab the ortho ax itself
                rSO = 1;
                for orthostep = setdiff(stepRangeOrtho, 0)

%                     if orthostep == 0
%                         featMat{cellIndex, realStepOrtho}(rSO, :) = -orthAx;
%                         featMat{cellIndex, realStepOrtho}(rSO+1, :) = orthAx;
%                         rSO = rSO + 2;
%                     else      
                        featVec = orthAx.*(std_params*orthostep); 
                        if compress_expand
                            featMat{cellIndex, realStepOrtho}(rSO, :) = featVec*V_inv + mu; 
                        else
                            featMat{cellIndex, realStepOrtho}(rSO, :) = featVec; 
                        end
                                               
                        rSO = rSO + 1;
%                     end
                end
                
            end
            
            % cell ID
            featMat{cellIndex, numcols+1} = responses{cellIndex, 3};
            cellIndex = cellIndex + 1;
        end
    end
end



% save([diskPath filesep taskPath filesep 'FeatureVectorsAlongAxis_NoNorm_ParamNorm_Scale5_' method type '_allITCells'], 'featMat');
if indiv_sess
    outPath = [diskPath filesep taskPath filesep 'predictedFeatures' filesep sess filesep image_spread];
else
    outPath = [diskPath filesep taskPath filesep 'predictedFeatures'];
end

if compress_expand
    outPath = [outPath filesep 'compress_expand'];
end

if ~exist(outPath, 'dir')
    mkdir(outPath)
end

%     if addNoise
%         save([outPath filesep ['FVAA_UnitVec_[1to' num2str(n_steps+1) ']_stepsize' num2str(scale) '_' method type '_' group '_' sess '_' image_spread '_' 'NoiseLevel5%']], 'featMat');
%     else
%         save([outPath filesep ['FVAA_UnitVec_[1to' num2str(n_steps+1) ']_stepsize' num2str(scale) '_' method type '_' group sess '_' image_spread '_' ]], 'featMat');
%     end

% Emily's python file reads each cell matrix separately
% cell array's in python
for cellIndex = 1:size(featMat, 1)
    if ~isempty(featMat{cellIndex, 1})
        for st = 1:size(featMat, 2)-1
        
            currMat = featMat{cellIndex, st};
            
            if strcmp(image_spread, 'PrefOnly')
                suffix = ['_' sprintf('%03d', st) '_(' num2str(stepRange(1)) 'to' num2str(stepRange(end)) ',0)_'];
            elseif strcmp(image_spread, 'PrefandOrtho')
                if st == 1
                    suffix = ['_' sprintf('%03d', st) '_(' num2str(stepRange(1)) 'to' num2str(stepRange(end)) ',0)_']; % pref axis (x axis)
                elseif st == 2
                    suffix = ['_' sprintf('%03d', st) '_(0,' num2str(stepRange(1)) 'to' num2str(stepRange(end)) ')_']; % orth axis (y axis)
                end
            elseif strcmp(image_spread, 'Grid')
                if st == size(featMat, 2)-1
                    suffix = ['_' sprintf('%03d', st) '_(0,' num2str(stepRange(1)) 'to' num2str(stepRange(end)) ')_'];
                else
                    suffix = ['_' sprintf('%03d', st) '_(' num2str(stepRange(1)) 'to' num2str(stepRange(end)) ',' num2str(stepRangeOrtho(st)) ')_'];
                end
            end
            
            save([outPath filesep ['cellid' num2str(strctCells(cellIndex).Name) '_cell' num2str(cellIndex) suffix 'stepsize' num2str(scale) '_' sess group ]], 'currMat');
            
        end
    end
end


%% Sanity check
% if ~strcmp(image_spread, 'PrefOnly')
%     
%     for i = 1:length(stepRange)       
%         u = featMat{2, end-1}(i, :);% orth ax
%         v = featMat{2, 1}(i, :); % pref ax
%         CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
%         ThetaInDegrees(i) = real(acosd(CosTheta));
%     end
%     
% end

%% helper function
function param = param_normalize_per_dim(param, amp_dim, NIMAGE)
%% normalize each dimension separately
%% Liang does this only - May2021
param = param./repmat(amp_dim, [NIMAGE 1]);
end



%% Useful code

% %% Sanity check - morph between 2 images
%
% % no need to compute axes, just take 2 imgs and produce points in between
%
% im1 = 49;
% im2 = 210;
%
% im1_p = params(im1, :);
% im2_p = params(im2, :);
%
% n_morph_steps = 25;
%
% m_stepsize = (im1_p - im2_p)/(n_morph_steps-1);
% s_beyond = 100;
% morphMat = [];
% ctr = 1;
% for m_s = -s_beyond:1:n_morph_steps+(s_beyond-1)
%
%     morphMat(ctr, :) = im2_p+((m_s)*m_stepsize);
%     ctr = ctr + 1;
% end
%
%
% save([diskPath filesep taskPath filesep 'predictedFeatures' filesep ['MorphSanityCheck_Ims' num2str(im1) '_' num2str(im2) '_' num2str(s_beyond) 'Beyond']], 'morphMat')
%
%
% %% How many dimensions is a cell tuned to?
%
% cellDim_tuning = [];
% alpha = 0.01;
% c = [];
% p = [];
%
% for cellIndex = 1:length(strctCells)
%     if cellIndex == 21
%         keyboard
%     end
%
%     fr1 = responses{cellIndex, 1}; % example cell 8276 object screening
%
%     for i = 1:size(params,2)
%         [cc pp] = corrcoef(params(:,i),fr1);
%         c(i) = cc(1,2); p(i) = pp(1,2);
%
%
%     end
%     ids = p < alpha;
%
%     cellDim_tuning(cellIndex, 1) = sum(double(ids));
%     %     cellDim_tuning(cellIndex, 1) = sum(double(ids));
%
% end
%
% f = figure;
% scatter(1:length(strctCells), cellDim_tuning, 'k', 'filled')
% title({['Number of dimensions a cell is tuned to: ' num2str(size(params, 2)) 'D space'], ['p < ' num2str(alpha)]})
% xlabel('Cell no')
% ylabel('No of dimensions')
% set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% filename = [diskPath filesep taskPath filesep 'predictedFeatures' filesep ['CellTuning_Dimensions_' num2str(size(params, 2)) 'D_Space']];
% print(f, filename, '-dpng', '-r0')
%
%
% %% compare the axes made by compression+expansion and the 4096D vector
%
%
% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_MatlabMean_' layermat '_' stimDir '.mat']]);
%
% method = 'sta';
% % method = 'linear_regression';
%
% for i = 1:2
%
%     switch i
%         case 1
%             params = feat;
%         case 2
%             [coeff, score, ~, ~, ~, mu] = pca(feat);%, 'Centered', false); % score is much higher range
%             params = score(:, 1:50);
%             V_inv = pinv(coeff(:, 1:50));
%             %             V_inv = coeff(:, 1:50)';
%
%     end
%
%     for cellIndex = 1:length(strctCells)
%         if ~isempty(responses{cellIndex, 2})
%
%             [sta, ~] = Utilities.ObjectSpace.analysis_STA(responses{cellIndex, 1}, params, method);
%
%             std_params = std(params); % std dev of each dimension
%
%             std_params(std_params == 0) = 1;
%
%             if strcmp(method, 'linear_regression')
%                 if strcmp(options.screenType, 'Face')
%
%                     % shape-appearance dimensions are flipped between mine and Stevens code
%                     % why is this only for lin-reg and not general sta?
%                     sta = fliplr(sta)./norm(sta); % normalize length
%
%                 else
%                     if i == 1
%                         sta = sta./norm(sta);
%                     elseif i == 2
%                         sta = sta*V_inv + mu;
%                         sta = sta./norm(sta);
%                     end
%                 end
%                 %             sta = sta*V_inv + mu;
%                 all_sta_linReg(cellIndex, :, i) = sta;
%
%             elseif strcmp(method, 'sta')
%                 if i == 1
%                     sta = sta./norm(sta);
%                 elseif i == 2
%                     sta = sta*V_inv + mu;
%                     sta = sta./norm(sta);
%                 end
%                 all_sta(cellIndex, :, i) = sta;
%
%             end
%         end
%     end
% end
%
%
%
% %%
% % all_sta = all_sta_linReg;
% t1 = all_sta(:, :, 1);
% t2 = all_sta(:, :, 2);
% ccc = diag(corr(t1', t2'));
%
%
% for j = 1:121
%     csim2(j, 1) = 1-pdist([t1(j, :); t2(j, :)], 'cosine');
% end
% csim = squareform(1-pdist([t1; t2], 'cosine'));
% csim = csim(length(strctCells)+1:end, 1:length(strctCells));
%
%
% f = figure;
% hold on
% set(gca, 'FontSize', 14, 'FontWeight', 'bold')
% scatter(1:length(strctCells), ccc, 'b', 'filled')
% scatter(1:length(strctCells), csim2, 'r', 'filled')
% title('Axis comparison 50D space vs 4096D space')
% xlabel('Cell no')
% lgnd = legend({'Corr', 'Cosine Sim'});
% filename = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'Axis_comparison_C&D_vs_OG'];
% print(f, filename, '-dpng', '-r0')


