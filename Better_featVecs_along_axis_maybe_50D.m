% /media/vwadia/T7/SUAnalysis/Localiser_Task/P85CS


% Better script for sampling images 50D


setDiskPaths
indiv_sess = 1;
compress_expand = 1;
splitByCat = 0;
normSpace = false;    


% taskPath = 'Recall_Task';
% sessPath = [diskPath filesep taskPath filesep 'P84CS' filesep 'RecallScreening_Session_2_20230408']; patID = 'P84CS';

taskPath = 'Object_Screening';
sessPath = [diskPath filesep taskPath filesep 'P85CS' filesep 'ClosedLoopScreening_Session_1_20230419']; patID = 'P85CS';

% taskPath = 'Recall_Task';
% sessPath = [diskPath filesep taskPath filesep 'P85CS' filesep 'RecallScreening_Session_1_20230424']; patID = 'P85CS';



load([sessPath filesep 'PsthandResponses']);
load([sessPath filesep 'strctCells']);
sess = strctCells(1).SessionID;
% load in params - adjust this for screen type

options.screenType = 'Object';
stimDir = '500Stimuli';
layermat = 'fc6';
options.ind_train = [1:500]';

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
python = true;

%% figure out max proj value

image_spread = 'PrefandOrtho';
% image_spread = 'Grid';
group = '';

if strcmp(patID, 'P84CS') && strcmp(image_spread, 'Grid')
    startCell = 4;
    lastCell = 4;
elseif strcmp(patID, 'P85CS') 
    if strcmp(taskPath, 'Object_Screening')
        startCell = 11;
        lastCell = 11;
    elseif strcmp(taskPath, 'Recall_Task')
        if strcmp(image_spread, 'Grid')
            startCell = 6;
            lastCell = 6;
        else
            startCell = 1;
            lastCell = 6;
        end
    end
else
    startCell = 1;
    lastCell = length(strctCells);
end

for cellIndex = startCell:lastCell
    
    
fr = responses{cellIndex, 1};

para = params(options.ind_train,:);
amp_dim = sqrt(sum(para.^2)); % finding the norm of each dimension 1xndim
amp_dim(amp_dim == 0) = 1;
if strcmp(options.screenType, 'Object')
    para = param_normalize_per_dim(para, amp_dim, length(options.ind_train));
elseif strcmp(options.screenType, 'Face')
    para = param_normalize(para, amp_dim, ndim1);
end
ndim = size(para, 2);
 
if isfield(options, 'sta') && ~isempty(options.sta)
    sta = options.sta;
else
    sta=fr'*para;
end

if normSpace
    value_sta_prj = (sta/norm(sta))*para';
else
    value_sta_prj = (sta/norm(sta))*params';
end

% scatter plot STA vs max orth STA
para_sub_sta = zeros(size(para));
for k=1:size(para,1);
    param_sta_prj = sta*(para(k,:)*sta')/(sta*sta'); % vector of params pojected onto STA
    para_sub_sta(k,:) = para(k,:) - param_sta_prj; % subtract STA component from param
end


if isfield(options, 'orthAx') && ~isempty(options.orthAx)
    COEFF(:, 1) = options.orthAx;
else
    % PCA
    COEFF = pca(para_sub_sta);
end

pc1 = para_sub_sta * COEFF(:,1);

% dot_color = zeros(size(fr,1),3); 
% dot_color(:,1) = ((fr-min_fr)/(max_fr-min_fr));
% dot_color(:,3) = 1- dot_color(:,1);

% h2 = subplot(3, 3, [5 6 8 9]);


[sorted_fr, reorder_ind] = sort(fr);
% [sorted_fr, reorder_ind] = sort(value_sta_prj); % is the same bc of how axis is computed
x = value_sta_prj(reorder_ind); 
y = pc1(reorder_ind);

maxProjPref = max(x);
minProjPref = min(x);

maxProjOrtho = max(y);
minProjOrtho = min(y);

%% compute step number required

if normSpace
    std_params = std(para);
else
    std_params = std(params);
end
sta = sta./norm(sta); % normalize length
all_sta(cellIndex, :) = sta;

stepLim = maxProjPref/ norm(sta.*(std_params));
edge = ceil(stepLim)+1;
numSteps = 10; % arbitrary choice
scale = (edge*2)/numSteps; 

stepRange = setdiff([-edge:scale:edge], 0);
stepRangeOrtho = stepRange;
orthAx = COEFF(:,1)';

%% sample along axis


if strcmp(image_spread, 'PrefandOrtho')
    
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
end

if indiv_sess
    outPath = [diskPath filesep taskPath filesep 'predictedFeatures' filesep sess filesep image_spread];
else
    outPath = [diskPath filesep taskPath filesep 'predictedFeatures'];
end

if compress_expand
    outPath = [outPath filesep 'compress_expand'];
end

if normSpace 
    outPath = [outPath filesep 'withNorm'];
end

if ~exist(outPath, 'dir')
    mkdir(outPath)
end

% save
for st = 1:size(featMat, 2)
    currMat = featMat{cellIndex, st};
    
    if strcmp(image_spread, 'PrefOnly')
        suffix = ['_' sprintf('%03d', st) '_(' num2str(-edge) 'to' num2str(edge) ',0)_'];
    elseif strcmp(image_spread, 'PrefandOrtho')
        if st == 1
            suffix = ['_' sprintf('%03d', st) '_(' num2str(-edge) 'to' num2str(edge) ',0)_']; % pref axis (x axis)
        elseif st == 2
            suffix = ['_' sprintf('%03d', st) '_(0,' num2str(-edge) 'to' num2str(edge) ')_']; % orth axis (y axis)
        end
    elseif strcmp(image_spread, 'Grid')
%         if st == size(featMat, 2)-1
%             suffix = ['_' sprintf('%03d', st) '_(0,' num2str(-edge) 'to' num2str(edge) ')_'];
%         else
            suffix = ['_' sprintf('%03d', st) '_(' num2str(-edge) 'to' num2str(edge) ',' num2str(stepRangeOrtho(st)) ')_'];
%         end
    end
    
    save([outPath filesep ['cellid' num2str(strctCells(cellIndex).Name) '_cell' num2str(cellIndex) suffix 'stepsize' num2str(scale) '_' sess group '.mat']], 'currMat');
end

end


%% saving for python
% for cellIndex = 1:size(featMat, 1)
%     if ~isempty(featMat{cellIndex, 1})
%         for st = 1:size(featMat, 2)-1
        
           
            
%         end
%     end
% end

%% HELPER
function param = param_normalize_per_dim(param, amp_dim, NIMAGE)
%% normalize each dimension separately 
%% Liang does this only - May2021
param = param./repmat(amp_dim, [NIMAGE 1]);
end
