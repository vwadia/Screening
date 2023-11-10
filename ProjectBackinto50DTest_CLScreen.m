

% projecting 4096D feature vectors of images used in CLScreen back into
% 50D and using in STA function to see spread.

% vwadia Jan2023



setDiskPaths


taskPath = 'Object_Screening';
m_sessPath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030']; patID = 'P81CS';
a_sessPath = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopReScreen_Session_1_20221030'];

% taskPath = 'Object_Screening';
% m_sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115'];patID = 'P82CS';
% a_sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'ClosedLoopReScreen_Session_1_20230115'];

% did not actually end up running this
% taskPath = 'Recall_Task';
% m_sessPath = [diskPath filesep 'Recall_Task' filesep 'P84CS' filesep 'RecallScreening_Session_2_20230408'];patID = 'P84CS';
% a_sessPath = [diskPath filesep 'Recall_Task' filesep 'P84CS' filesep 'RecallScreening_Session_2_20230408'];

% taskPath = 'Object_Screening';
% m_sessPath = [diskPath filesep taskPath filesep 'P85CS' filesep 'ClosedLoopScreening_Session_1_20230419'];patID = 'P85CS';
% a_sessPath = [diskPath filesep taskPath filesep 'P85CS' filesep 'ClosedLoopReScreen_Session_1_20230419'];

% taskPath = 'Recall_Task';
% m_sessPath = [diskPath filesep taskPath filesep 'P85CS' filesep 'RecallScreening_Session_1_20230424'];patID = 'P85CS';
% a_sessPath = [diskPath filesep taskPath filesep 'P85CS' filesep 'ReScreenRecall_Session_1_20230424'];


stimPath = [a_sessPath filesep 'genStimOnly'];

if ~strcmp(patID, 'P84CS')
    load([a_sessPath filesep 'SynthPsthandResponses'])
end
a_strctCells = load([a_sessPath filesep 'strctCells']);

load([a_sessPath filesep 'PsthandResponses'])

m_strctCells = load([m_sessPath filesep 'strctCells']);


%%
clearvars options



if strcmp(patID, 'P81CS')
    
    n_steps = 10;
    n_steps_ortho = 10; % set per session
    
    greyImsOnly = 0; % only black and white images
    
    corrThresh = 0.9;
    morn_cells = [1100; 1546; 1877; 1897; 2156];
    aft_matches = [1525; 1943; 1928; 2300; 1935];
    aft_idx = [1; 2; 3; 4; 5];
    legNames = {'1525', '1943', '1928', '2300', '1935'};
    
    GridOnly = 0; % by default
    dimRed = 1; % if yes do dimensionality reduction and plot in 50D space otherwise do in 4096D space
    useBigAxes = 0; % compute axes in 4096 space and convert them to 50D space (instead of recomputing in 50D)
    
elseif strcmp(patID, 'P82CS')
    
    corrThresh = 0.9;
    morn_cells = [nan; nan; 904]; % placeholding so cellIndex is correct
    aft_matches = [nan; nan; 2360];
    aft_idx = [1];
    
    n_steps = 5;
    n_steps_ortho = 5; % set per session
    
    greyImsOnly = 0; % only black and white images
    GridOnly = 0;
    
    % to recreate thesis figure = dimRed = 0 and useBigAxes = 1 w/ morning axis
    dimRed = 1; % if yes do dimensionality reduction and plot in 50D space otherwise do in 4096D space
    useBigAxes = 0; % compute axes in 4096 space and convert them to 50D space (instead of recomputing in 50D)
    
    
elseif strcmp(patID, 'P84CS')
    
    corrThresh = 0.9;
    morn_cells = nan(length(m_strctCells.strctCells), 1);
    morn_cells(2, 1) = 1482; morn_cells(4, 1) = 2040;
    aft_matches = morn_cells;
    aft_idx = 1;
    
    n_steps = 5;
    n_steps_ortho = 5; % set per session
    
    greyImsOnly = 0; % only black and white images
    GridOnly = 0;
    
    dimRed = 1; % if yes do dimensionality reduction and plot in 50D space otherwise do in 4096D space
    useBigAxes = 0; % compute axes in 4096 space and convert them to 50D space (instead of recomputing in 50D)

elseif strcmp(patID, 'P85CS')
    
    if strcmp(taskPath, 'Object_Screening')
        corrThresh = 0.9;
        morn_cells = nan(length(m_strctCells.strctCells), 1);
        morn_cells(7, 1) = 1355;
        morn_cells(11, 1) = 1419;
        aft_matches = morn_cells;
        aft_matches(7, 1) = 1580;
        aft_matches(11, 1) = 1539; % not actually a match
        aft_idx = 1;
        
        n_steps = 5;
        n_steps_ortho = 5; % set per session
        
        greyImsOnly = 0; % only black and white images
        GridOnly = 0;
        
        dimRed = 1; % if yes do dimensionality reduction and plot in 50D space otherwise do in 4096D space
        useBigAxes = 0; % compute axes in 4096 space and convert them to 50D space (instead of recomputing in 50D)

    elseif strcmp(taskPath, 'Recall_Task')
        corrThresh = 0.9;
        morn_cells = cat(1, m_strctCells.strctCells(:).Name);
        morn_cells(3, 1) = nan; morn_cells(4, 1) = nan;
        aft_matches = [8320; 8302; nan; nan; 4305; 1189]; % original
        
        aft_idx = 1;
        
        n_steps = 5;
        n_steps_ortho = 5; % set per session
        
        greyImsOnly = 0; % only black and white images
        GridOnly = 0;
        
        dimRed = 1; % if yes do dimensionality reduction and plot in 50D space otherwise do in 4096D space
        useBigAxes = 0; % compute axes in 4096 space and convert them to 50D space (instead of recomputing in 50D)
        
    end
    
end

% xlim
fullImDir = Utilities.readInFiles(stimPath);
fullImDirCell = struct2cell(fullImDir)';
if strcmp(patID, 'P85CS') && strcmp(taskPath, 'Recall_Task')
    [~, natIdx] = natsortfiles({fullImDir.name});
    fullImDirCell = fullImDirCell(natIdx, 1);
    fullImDir = fullImDir(natIdx, 1);
end

if strcmp(patID, 'P81CS')
    startCell = 1;
    endCell = length(morn_cells);
elseif strcmp(patID, 'P82CS')
    startCell = 3;
    endCell = length(morn_cells);
elseif strcmp(patID, 'P84CS')
    startCell = 2;
    endCell = 2;
elseif strcmp(patID, 'P85CS')
    if strcmp(taskPath, 'Object_Screening')
        startCell = 7;
        endCell = 11;
    elseif strcmp(taskPath, 'Recall_Task')
        startCell = 1;
        endCell = 6;
    end
end

for cellIndex = startCell:endCell
    if ~isnan(morn_cells(cellIndex))
        mornCell = morn_cells(cellIndex);
        
        % images corresponding to that cell
        cellIms = cellfun(@(x) strcmp(x(7:7+numel(num2str(mornCell))-1), num2str(mornCell)), fullImDirCell(:, 1), 'UniformOutput', false);
        idx = structfind(a_strctCells.strctCells, 'Name', aft_matches(cellIndex));
        if ~strcmp(patID, 'P84CS')
            resp = synthResponses{idx, 1}(cell2mat(cellIms), 1); % responses only to desired images
        end
        imNames = fullImDirCell(cell2mat(cellIms), 1);
        
        if GridOnly
            dashpos = strfind(imNames{1}, '_');
            cellGridIms = cellfun(@(x) strcmp(x(dashpos(1)+1:dashpos(1)+4), 'Grid'), imNames, 'UniformOutput', false);
            imNames = imNames(cell2mat(cellGridIms));
            % only responses to grid images (for now)
            resp = resp(cell2mat(cellGridIms));
        else
            if strcmp(patID, 'P82CS')
                imNames = imNames(1:40);
                resp = resp(1:40);
            elseif strcmp(patID, 'P85CS')
                if strcmp(taskPath, 'Object_Screening')
                    imNames = imNames(1:20);
                    resp = resp(1:20);
                elseif strcmp(taskPath, 'Recall_Task')
                    if isequal(mornCell, 2700) % only the pref and ortho from the grid cell
                        imNames = imNames(1:20);
                        resp = resp(1:20);
                    end
                end
            end           
        end
        
        if strcmp(patID, 'P82CS')
            % b&W images came first
            if greyImsOnly
                imNames = imNames(1:length(imNames)/2);
                resp = resp(1:length(resp)/2);
            else
                imNames = imNames((length(imNames)/2)+1:end);
                resp = resp((length(resp)/2)+1:end);
            end
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
        
        if strcmp(patID, 'P81CS')
            
            % bootleg way for now ---------------------------------------------
            featPath = [diskPath filesep taskPath filesep 'predictedFeatures'];
            load([featPath filesep 'FVAA_UnitVec_[1to10]_stepsize1_sta_P81CS_AM']); % loads featMat n_cells x 2 cell array with matrices and cell ids
            
            c_idx = cellfun(@(x) isequal(morn_cells(cellIndex), x), featMat(:, 2), 'UniformOutput', false);
            c_idx = cell2mat(c_idx);
            
            currMat = featMat{c_idx};
            prefAxIms = currMat;
            options.famNorm = false;
            % -----------------------------------------------------------------
            
            % proper way - re-save features so it works like this March2023
            % vwadia
            % --------------------------------------------------------------
            %         featPath = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'P81CS_AM'];
            
            % load preferred axis
            %         load([featPath filesep 'PrefandOrtho' filesep 'cellid' num2str(morn_cells(cellIndex)) '_cell6_001_(-10to10,0)_stepsize2_P81CS_AM']); type = 'PrefandOrtho';
            
            % load orthogonal axis - Where did this come from? I never ran the
            % ortho axis in P81...
            %         load([featPath filesep 'PrefandOrtho' filesep 'cellid' num2str(morn_cells(cellIndex)) '_cell6_001_(0,-10to10)_stepsize2_P81CS_AM']);
            % -----------------------------------------------------------
        elseif strcmp(patID, 'P82CS')
            featPath = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'P82CS_CL_1'];
            
            % load preferred axis
            if GridOnly
                options.famNorm = false;
                load([featPath filesep 'Grid' filesep 'cellid904_cell3_006_(-10to10,0)_stepsize2_P82CS_CL_1']); type = 'Grid';
            else
                options.famNorm = false;
                load([featPath filesep 'PrefandOrtho' filesep 'cellid904_cell3_001_(-50to50,0)_stepsize10_P82CS_CL_1']); type = 'PrefandOrtho';
            end
            prefAxIms = currMat;
            
            % load orthogonal axis
            if GridOnly
                options.unfamNorm = false;
                load([featPath filesep 'Grid' filesep 'cellid904_cell3_012_(0,-10to10)_stepsize2_P82CS_CL_1']);
            else
                options.unfamNorm = false;
                load([featPath filesep 'PrefandOrtho' filesep 'cellid904_cell3_002_(0,-50to50)_stepsize10_P82CS_CL_1']);
            end
            orthAxIms = currMat;
            
        elseif strcmp(patID, 'P84CS')
            featPath = [diskPath filesep 'Recall_Task' filesep 'predictedFeatures' filesep 'P84CS_RecScreen_2'];
            
            % load preferred axis
            if GridOnly
                keyboard % not relevant for this patient
            else
                options.famNorm = false;
                load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid1482_cell2_001_(-6to6,0)_stepsize1.2_P84CS_RecScreen_2.mat']); type = 'PrefandOrtho';
            end
            prefAxIms = currMat;
            
            % load orthogonal axis
            if GridOnly
                keyboard % not relevant for this patient
            else
                options.unfamNorm = false;
                load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid1482_cell2_002_(0,-6to6)_stepsize1.2_P84CS_RecScreen_2.mat']); type = 'PrefandOrtho';
            end
            orthAxIms = currMat;
            
        elseif strcmp(patID, 'P85CS')
            if strcmp(taskPath, 'Object_Screening')
                featPath = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'P85CS_1'];
                
                % load preferred axis
                if GridOnly
                    keyboard
                else
                    options.famNorm = false;
                    if cellIndex == 11
                        load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid1419_cell11_001_(-6to6,0)_stepsize1.2_P85CS_1.mat']); type = 'PrefandOrtho';
                    elseif cellIndex == 7
                        load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid1355_cell7_001_(-4to4,0)_stepsize0.8_P85CS_1.mat']); type = 'PrefandOrtho';
                    end
                end
                prefAxIms = currMat;
                
                % load orthogonal axis
                if GridOnly
                    keyboard
                else
                    options.unfamNorm = false;
                    if cellIndex == 11
                        load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid1419_cell11_002_(0,-6to6)_stepsize1.2_P85CS_1.mat']);
                    elseif cellIndex == 7
                        load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep 'cellid1355_cell7_002_(0,-4to4)_stepsize0.8_P85CS_1.mat']);
                    end
                end
                orthAxIms = currMat;
                
            elseif strcmp(taskPath, 'Recall_Task')
                featPath = [diskPath filesep taskPath filesep 'predictedFeatures' filesep 'P85CS_RecScreen_1'];
                
                % load preferred axis
                if GridOnly
                    keyboard
                else
                    options.famNorm = false;
                    edges = [5; 6; 5; 5; 5; 5];
                    ssize = [1; 1.2; 1; 1; 1; 1];
                    if ~isnan(morn_cells(cellIndex))
                        featname = ['cellid' num2str(morn_cells(cellIndex)) '_cell' num2str(cellIndex)...
                            '_001_(-' num2str(edges(cellIndex)) 'to' num2str(edges(cellIndex)) ',0)_stepsize' num2str(ssize(cellIndex)) '_P85CS_RecScreen_1.mat'];
                        disp(featname);
                        load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep featname]); type = 'PrefandOrtho';
                    end
                end
                prefAxIms = currMat;
                
                % load orthogonal axis
                if GridOnly
                    keyboard
                else
                    options.unfamNorm = false;
                    edges = [5; 6; 5; 5; 5; 5];
                    ssize = [1; 1.2; 1; 1; 1; 1];
                    if ~isnan(morn_cells(cellIndex))
                        featname = ['cellid' num2str(morn_cells(cellIndex)) '_cell' num2str(cellIndex)...
                            '_002_(0,-' num2str(edges(cellIndex)) 'to' num2str(edges(cellIndex)) ')_stepsize' num2str(ssize(cellIndex)) '_P85CS_RecScreen_1.mat'];
                        disp(featname);
                        load([featPath filesep 'PrefandOrtho' filesep 'compress_expand' filesep featname]); type = 'PrefandOrtho';
                    end
                end
                orthAxIms = currMat;
            end
        end
        
        if ~strcmp(patID, 'P81CS')
            axIms = [prefAxIms; orthAxIms];
        else
            axIms = [prefAxIms];
        end
        
        if dimRed
            proj_into_500 = (axIms - repmat(mu, [size(axIms, 1) 1]))*coeff; % coeff = PCs of 500 object space
            proj = proj_into_500(:, 1:50); % now the images are projected into the space built by my 500 ims
        else
            proj = axIms;
        end
        
        % collect responses to those ims only
        if ~strcmp(patID, 'P84CS')
            rr = resp;
            rr_ortho = resp(find(x == 0));
            rr_pref = resp(find(y == 0));
        end
        %
        % compute axis
        load([m_sessPath filesep 'PsthandResponses']); mornAxis = true;
%         load([a_sessPath filesep 'PsthandResponses']); oldcellIndex = cellIndex; cellIndex = idx;  mornAxis = false;
        
        if useBigAxes && dimRed
            
            [sta, ~] = Utilities.ObjectSpace.analysis_STA(responses{cellIndex, 1}, feat, 'sta');
            
            staProj500 = (sta - mu)*coeff;
            %         staProj500 = (sta - mean(sta))*coeff;
            %         staProj500 = (sta - mean(axIms, 1))*coeff;
            
            options.sta = staProj500(1, 1:size(params, 2));
            sta = options.sta;
            % compute orthax in 4096D
%             para = feat(options.ind_train,:);
            para = params(options.ind_train,:);
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
            
%             orthAxProj500 = (orthAx - mu)*coeff;
            %         orthAxProj500 = (orthAx - mean(orthAx))*coeff;
            %         orthAxProj500 = (orthAx - mean(axIms, 1))*coeff;
            
%             options.orthAx = orthAxProj500(1, 1:size(params, 2))';
            options.orthAx = orthAx; % uncomment this vwadia 9/14
            
            
        else
            [sta, ~] = Utilities.ObjectSpace.analysis_STA(responses{cellIndex, 1}, params, 'sta');
            options.sta = sta;
        end
        
        params = [params; proj];
        options.fam = 1;
        num_fam = n_steps*2; % both directions
        options.fam_stim_ind = 500+1:500+num_fam;
        options.fam_para_ind = 500+1:500+num_fam;
        
        if exist('orthAxIms', 'var')
            num_unfam = n_steps_ortho*2; % both directions
            options.unfam = 1;
            options.unfam_stim_ind = 500+num_fam+1:500+num_fam+num_unfam;
            options.unfam_para_ind = 500+num_fam+1:500+num_fam+num_unfam;
        end
        if ~strcmp(patID, 'P84CS')
            allResp = [responses{cellIndex, 1}; rr_pref; rr_ortho];
        else
            allResp = [responses{cellIndex, 1}; ones(n_steps*2, 1); ones(n_steps_ortho*2, 1)];
        end
        
        hfig = Utilities.ObjectSpace.STA_figure_original(allResp, params, options);
        
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontWeight','Bold', 'LineWidth', 1.2);

        
        ndim = size(params, 2);
        
        if greyImsOnly
            filename = [a_sessPath filesep 'SynthIms_STAPlot_' num2str(ndim) 'D'];
        else
            filename = [a_sessPath filesep 'SynthImsColour_STAPlot_' num2str(ndim) 'D'];
        end
        if mornAxis
            filename = [filename '_mornAxis_' num2str(morn_cells(cellIndex)) '_' num2str(aft_matches(cellIndex))];
        else
            filename = [filename '_aftAxis_' num2str(morn_cells(oldcellIndex)) '_' num2str(aft_matches(oldcellIndex))];            
        end
%             keyboard
%         print(hfig, filename, '-dpng', '-r0');
%         close all
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


