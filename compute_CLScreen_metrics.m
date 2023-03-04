% script to compute other metrics related to CLScreen

% 1 - The dynamic range of each neuron
%     Per cell:
%     compute axis - project all images and choose max and min values
%     how many steps does it take in each direction to get there?
     
% 2 - comparing the firng rate of STA image to baseline for all cells - does the image being colored affect the overall FR?

% 3 - are the synth images off manifold wrt to the stim set? See Liang's paper Fig 2a Fam vs Unfam

%% set paths and load data

setDiskPaths

% a_basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030'];
% load([a_basePath filesep 'SynthPsthandResponses.mat'])
% load([a_basePath filesep 'strctCells.mat'])
% 
% m_basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030'];
% morn = load([m_basePath filesep 'PsthandResponses']);
% aft = load([a_basePath filesep 'PsthandResponses']);
compress_expand = 1;

taskPath = 'Object_Screening';
sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'FingerprintScreening_Session_1_20230111'];

load([sessPath filesep 'PsthandResponses']);
load([sessPath filesep 'strctCells']);
sess = strctCells(1).SessionID;

imageIDs = [1:500]';
options.screenType = 'Object';
stimDir = '500Stimuli';
layermat = 'fc6';


load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_MatlabMean_' layermat '_' stimDir '.mat']]); % used for CLScreen
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


%% compute dynamic range (as check/confirmation for the afternoon cells)

% range = (max - min projection value)/std_dev of dimensions  

% compute axis for each
method = 'sta';
% method = 'linear_regression';

% normalize parameters
para = params;
amp_dim = sqrt(sum(para.^2)); % finding the norm of each dimension 1xndim
amp_dim(amp_dim == 0) = 1; % this makes no difference when it is 50 dim

if strcmp(options.screenType, 'Object')
    para = param_normalize_per_dim(para, amp_dim, length(imageIDs));
elseif strcmp(options.screenType, 'Face')
    para = param_normalize(para, amp_dim, ndim1);
end

std_params = std(params);
std_para = std(para);
ndim = size(para, 2);



for cellIndex = 1:length(strctCells)
    
    if ~isempty(responses{cellIndex, 2})
        
        resp = responses{cellIndex, 1};
        fr = resp;
        % compute axis - feed in regular params as they are normalized
        % inside function
        [sta, ~] = Utilities.ObjectSpace.analysis_STA(resp, params, method);
        
        para_sub_sta = zeros(size(para));
        params_sub_sta = zeros(size(params));
        for k=1:size(para,1)
            param_sta_prj = sta*(para(k,:)*sta')/(sta*sta'); % vector of params pojected onto STA
            para_sub_sta(k,:) = para(k,:) - param_sta_prj; % subtract STA component from param          
        end
        
        % arrange ortho axes by how much variability they explain
        COEFF = pca(para_sub_sta);
        
        % should I use non-normalized parameters for this??
        value_ortho_prj = para_sub_sta * COEFF(:,1);

        value_sta_prj = (sta/norm(sta))*para';
                
       
        

   
    end
end

%% helpers
function param = param_normalize(param, amp_dim, ndim1)
%% normalize shape/appearance separately while keeping the relative amplitude within shape or appearance dimensions
%% stevens way - in the cell paper
ndim = size(param, 2);
% para = para./repmat(amp_dim, [NIMAGE 1]);

param(:,1:ndim1)=param(:,1:ndim1) / sqrt(sum(amp_dim(1:ndim1).^2)) / sqrt(2);
param(:, ndim1+1:ndim)=param(:, ndim1+1:ndim) / sqrt(sum(amp_dim(ndim1+1:ndim).^2)) / sqrt(2);
end

function param = param_normalize_per_dim(param, amp_dim, NIMAGE)
%% normalize each dimension separately 
%% Liang does this only - May2021
param = param./repmat(amp_dim, [NIMAGE 1]);
end

