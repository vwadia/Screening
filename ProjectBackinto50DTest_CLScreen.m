

% projecting 4096D feature vectors of images used in CLScreen back into
% 50D and using in STA function to9 see spread.

% vwadia Jan2023



setDiskPaths

taskPath = 'Object_Screening';


% m_basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030'];patID = 'P81CS';
% a_basePath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030'];

m_sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115'];patID = 'P82CS';
a_sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'ClosedLoopReScreen_Session_1_20230115'];


stimPath = [a_sessPath filesep 'genStimOnly'];

load([a_sessPath filesep 'SynthPsthandResponses'])
a_strctCells = load([a_sessPath filesep 'strctCells']);

load([a_sessPath filesep 'PsthandResponses'])


m_strctCells = load([m_sessPath filesep 'strctCells']);


%%

dimRed = 1; % if yes do dimensionality reduction and plot in 50D space otherwise do in 4096D space

useBigAxes = 1; % compute axes in 4096 space and convert them to 50D space (instead of recomputing in 50D)

if strcmp(patID, 'P82CS')
    
    corrThresh = 0.9;
    morn_cells = [nan; nan; 904]; % placeholding so cellIndex is correct
    aft_matches = [nan; nan; 2360]; 
    aft_idx = [1];

    n_steps = 5;
    n_steps_ortho = 5; % set per session
    scale = 2;
    stepRangeOrtho = [-n_steps_ortho:1:n_steps_ortho]*scale;
    stepRange = [-n_steps:1:n_steps]*scale;
    nbin = 10;

    greyImsOnly = 0; % only black and white images
    GridOnly = 0;

elseif strcmp(patID, 'P81CS')

    n_steps = 10;
    n_steps_ortho = 10; % set per session
    scale = 1;
    stepRangeOrtho = [-n_steps_ortho:1:n_steps_ortho]*scale;
    stepRange = [-n_steps:1:n_steps]*scale;
    nbin = 10;

    greyImsOnly = 0; % only black and white images

    corrThresh = 0.9;
    morn_cells = [1100; 1546; 1877; 1897; 2156];
    aft_matches = [1525; 1943; 1928; 2300; 1935];
    aft_idx = [1; 2; 3; 4; 5];
    legNames = {'1525', '1943', '1928', '2300', '1935'};

    GridOnly = 0; % by default

end


fullImDir = Utilities.readInFiles(stimPath);
fullImDirCell = struct2cell(fullImDir)';

for cellIndex = 3

    mornCell = morn_cells(cellIndex);

    % images corresponding to that cell
    cellIms = cellfun(@(x) strcmp(x(7:7+numel(num2str(mornCell))-1), num2str(mornCell)), fullImDirCell(:, 1), 'UniformOutput', false);
    idx = structfind(a_strctCells.strctCells, 'Name', aft_matches(cellIndex));
    resp = synthResponses{idx, 1}(cell2mat(cellIms), 1); % responses only to desired images
    imNames = fullImDirCell(cell2mat(cellIms), 1);

    

    if GridOnly
        cellGridIms = cellfun(@(x) strcmp(x(dashpos(1)+1:dashpos(1)+4), 'Grid'), imNames, 'UniformOutput', false);
        imNames = imNames(cell2mat(cellGridIms));
        % only responses to grid images (for now)
        resp = resp(cell2mat(cellGridIms));
    else
        imNames = imNames(1:40);
        resp = resp(1:40);
    end

    % b&W images came first
    if greyImsOnly
        imNames = imNames(1:length(imNames)/2);
        resp = resp(1:length(resp)/2);
    else
        imNames = imNames((length(imNames)/2)+1:end);
        resp = resp((length(resp)/2)+1:end);
    end

    for im = 1:length(imNames)

        brapos = strfind(imNames{im}, '(');
        ketpos = strfind(imNames{im}, ')');
        compos = strfind(imNames{im}, ',');

        x(im) = str2num(imNames{im}(brapos+1:compos-1));
        y(im) = str2num(imNames{im}(compos+1:ketpos-1));

    end


    stimDir = '500Stimuli';
    layermat = 'fc6';
    load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_MatlabMean_' layermat '_' stimDir '.mat']]);
    options.screenType = 'Object';
    options.ind_train = [1:500]';
    

    if dimRed 
        [coeff, score, ~, ~, ~, mu] = pca(feat);%, 'Centered', false); % score is much higher range
        params = score(:, 1:50);
        V_inv = pinv(coeff(:, 1:50));
    else
        params = feat;
    end

    featPath = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'P82CS_CL_1'];

    % load preferred axis
    % load([featPath filesep 'Grid' filesep 'cellid904_cell3_006_(-10to10,0)_stepsize2_P82CS_CL_1']); type = 'Grid';
%     load([featPath filesep 'PrefandOrtho' filesep 'cellid904_cell3_001_(-50to50,0)_stepsize10_P82CS_CL_1']); type = 'PrefandOrtho';
%     load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid904_cell3_001_(-10to10,0)_stepsize2_P82CS_CL_1']); type = 'PrefandOrtho';
    load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid904_cell3_001_(-5to5,0)_stepsize1_P82CS_CL_1']); type = 'PrefandOrtho';
    prefAxIms = currMat;

    % load orthogonal axis
    % load([featPath filesep 'Grid' filesep 'cellid904_cell3_012_(0,-10to10)_stepsize2_P82CS_CL_1']);
%     load([featPath filesep 'PrefandOrtho' filesep 'cellid904_cell3_002_(0,-50to50)_stepsize10_P82CS_CL_1']);
%     load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid904_cell3_002_(0,-10to10)_stepsize2_P82CS_CL_1']);
    load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid904_cell3_002_(0,-5to5)_stepsize1_P82CS_CL_1']);
    orthAxIms = currMat;

    axIms = [prefAxIms; orthAxIms];
    % axIms = [prefAxIms];
    % axIms = [orthAxIms];
    if dimRed
        proj_into_500 = (axIms - repmat(mu, [size(axIms, 1) 1]))*coeff; % coeff = PCs of 500 object space
        proj = proj_into_500(:, 1:50); % now the images are projected into the space built by my 500 ims
    else
        proj = axIms;
    end

    % collect responses to those ims only
    rr = resp;
    rr_ortho = resp(find(x == 0));
    rr_pref = resp(find(y == 0));

    % compute axis
    load([m_sessPath filesep 'PsthandResponses'])

    if useBigAxes && dimRed

        [sta, ~] = Utilities.ObjectSpace.analysis_STA(responses{cellIndex, 1}, feat, 'sta');
        
        staProj500 = (sta - mu)*coeff;
%         staProj500 = (sta - mean(sta))*coeff;
%         staProj500 = (sta - mean(axIms, 1))*coeff;
        
        options.sta = staProj500(1, 1:size(params, 2));

        % compute orthax in 4096D
        para = feat(options.ind_train,:);
        amp_dim = sqrt(sum(para.^2)); % finding the norm of each dimension 1xndim
        amp_dim(amp_dim == 0) = 1;
        % normalize features
        para_sub_sta = zeros(size(para));
        if strcmp(options.screenType, 'Object')
            para = param_normalize_per_dim(para, amp_dim, length(options.ind_train));
        elseif strcmp(options.screenType, 'Face')
            para = param_normalize(para, amp_dim, ndim1);
        end
        % subtract STA component
        for k=1:size(para,1);
            param_sta_prj = sta*(para(k,:)*sta')/(sta*sta'); % vector of params pojected onto STA
            para_sub_sta(k,:) = para(k,:) - param_sta_prj; % subtract STA component from param
        end
        % pick axis w/ most variability
        COEFF = pca(para_sub_sta);
        orthAx = COEFF(:, 1)';

        orthAxProj500 = (orthAx - mu)*coeff;
%         orthAxProj500 = (orthAx - mean(orthAx))*coeff;
%         orthAxProj500 = (orthAx - mean(axIms, 1))*coeff;

        options.orthAx = orthAxProj500(1, 1:size(params, 2))';



    else
        [sta, ~] = Utilities.ObjectSpace.analysis_STA(responses{cellIndex, 1}, params, 'sta');
        options.sta = sta;
    end
    
    params = [params; proj];
    options.fam = 1;
    options.fam_stim_ind = 501:510;
    options.fam_para_ind = 501:510;
    options.unfam = 1;
    options.unfam_stim_ind = 511:520;
    options.unfam_para_ind = 511:520;
    allResp = [responses{cellIndex, 1}; rr_pref; rr_ortho];
    % allResp = [responses{cellIndex, 1}; rr_pref];
    % allResp = [responses{cellIndex, 1}; rr_ortho];

    hfig = Utilities.ObjectSpace.STA_figure_original(allResp, params, options);
    
    ndim = size(params, 2);
   
    if greyImsOnly
        filename = [a_sessPath filesep 'SynthIms_STAPlot_' num2str(ndim) 'D'];
    else
        filename = [a_sessPath filesep 'SynthImsColour_STAPlot_' num2str(ndim) 'D'];
    end
%     print(hfig, filename, '-dpng', '-r0');
%     close all

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


