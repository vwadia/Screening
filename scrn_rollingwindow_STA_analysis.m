
%% 1. load in the data

% load in screeningData
% load in strctCells
% setTTLCodes

%% 2. set up rolling window and make matrices

firstWindow = -screeningData.timelimits(1)*1e3 + 50; % ms post stimON
lastWindow = -screeningData.timelimits(1)*1e3 + screeningData.timelimits(2)*1e3 - ceil(stimDur);
stepSize = 25;

% collect responses for each category in appropriate time window
for cellIndex = l(strctCells)
    ctr = 1;
    for window = firstWindow:stepSize:lastWindow
        windowLength = floor(stimDur);
        
        for i = l(screeningData.imageIDs)
            stimRas = screeningData.psth{cellIndex, 1}(find(screeningData.sortedOrder == screeningData.imageIDs(i)), window:window+windowLength);
            screeningData.rollingWindowSTA{cellIndex, 1}(i, ctr) = mean(mean(stimRas))*1e3;
        end
        screeningData.rollingWindowSTA{cellIndex, 2} = strctCells(cellIndex).Name;
        ctr = ctr+1;
        
    end
end
%%
if ~exist('params')
    if strcmp(taskStruct.subID, 'P73CS_Full')
        %     z_scored = 0;
        load([diskPath filesep 'ObjectSpace' filesep 'parameters_2k_synthetic_faces.mat']);
        params = params(1:667, :);
        options.screenType = 'Face';
    elseif strcmp(taskStruct.subID, 'P73CS_ParamObj') || strcmp(taskStruct.subID, '62')
        load([diskPath filesep 'ObjectSpace' filesep 'parameters_1593_objects.mat']); % will create score = 1593x50
        params = score;
        options.screenType = 'Object';
    else
        % get params matrix
         
        params = getDNActivations(pathStimuli, screeningData.imageIDs, 50, 'fc6');
        options.screenType = 'Object';
    end
end
pathOut = [basePath filesep 'STA_and_projections_rollingWindow_' num2str(stepSize) 'ms'];
if ~exist(pathOut)
    mkdir([pathOut]);
end
%% make STA plots
% for cellIndex = l(strctCells)
%     ctr = 1;
%     if ~isempty(screeningData.responses{cellIndex, 2})  
%         for window = firstWindow:stepSize:lastWindow
%             realWin = window + (screeningData.timelimits(1)*1e3);
%             options.ind_train = screeningData.imageIDs; % use all objects to calculate STA
%             [hfig, p] = STA_figure_original(screeningData.rollingWindowSTA{cellIndex, 1}(:,ctr), params, options);
%             sgtitle({['Cell number ' num2str(strctCells(cellIndex).Name)] 'STA and projections ', [strctCells(cellIndex).brainArea '\_Window\_' num2str(realWin)]});
%             
%             if p < 0.01
%                 newPathOut = [pathOut filesep 'significant_cells'];
%                 if ~exist(newPathOut)
%                     mkdir([newPathOut]);
%                 end
%                 if isfield(options, 'encoded_stim')
%                     print(hfig, [newPathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_orderedStim'], '-dpng', '-r0')
%                 else
%                     print(hfig, [newPathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name)], '-dpng', '-r0')
%                 end
%             else
%                 print(hfig, [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_' num2str(realWin)], '-dpng', '-r0')
%             end
%             
%             close all;
%             ctr = ctr+1;
%         end
%     end
% end

%% computing quality of rampiness - pearson/spearman OR just EV in linear model

% EV computation
for cellIndex = l(strctCells)
    ctr = 1;        
    disp(cellIndex);

    if ~isempty(screeningData.responses{cellIndex, 2}) % responsive cells only
        for window = firstWindow:stepSize:lastWindow
            realWin = window + (screeningData.timelimits(1)*1e3);           
            
            % observed responses
            resp = screeningData.rollingWindowSTA{cellIndex, 1}(:, ctr);
            
            % predicted responses
            [pred_resp, obs_resp] = Code.Utilities.computePredictedResponses(resp, params, screeningData.imageIDs, options.screenType);
%             screeningData.pred_resp{cellIndex, 1} = pred_resp;
            
            % calculate explained variance
            screeningData.evRollingWindow(cellIndex, ctr) = Code.Utilities.computeExplainedVariance(obs_resp, pred_resp);
            ctr = ctr+1;
        end
        [maxEVCell, pos_maxEVCell] = max(screeningData.evRollingWindow(cellIndex,:));
        windowRange = firstWindow+(screeningData.timelimits(1)*1e3):stepSize:lastWindow+(screeningData.timelimits(1)*1e3);
        yc_ev = find(windowRange == screeningData.responses{cellIndex, 2});
        if isequal(yc_ev, pos_maxEVCell)
            screeningData.rollingWindowSTA{cellIndex, 3} = 1;
        else
            screeningData.rollingWindowSTA{cellIndex, 3} = 0;
        end

        ispos = find(screeningData.evRollingWindow(cellIndex, :) > 0);
        if ~isempty(ispos)
            % plot relationship
            f = figure; clf;
            hold on
            plot(firstWindow+(screeningData.timelimits(1)*1e3):stepSize:lastWindow+(screeningData.timelimits(1)*1e3), screeningData.evRollingWindow(cellIndex, :));
            if ~isempty(yc_ev)
                plot(screeningData.responses{cellIndex, 2}, screeningData.evRollingWindow(cellIndex, yc_ev), '*','MarkerSize',10);
            end
            set(gca,'FontWeight', 'bold')
            yline(0, '--',  'FontWeight', 'bold');
            title({['EV vs Rolling Window'], [strctCells(cellIndex).brainArea '\_' num2str(strctCells(cellIndex).ChannelNumber) '\_' num2str(strctCells(cellIndex).Name)]});
            xlabel('WindowBegin post stimON', 'FontWeight', 'bold');
            ylabel('ev linear model', 'FontWeight', 'bold');
            print(f, [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_evWindow'], '-dpng', '-r0');
            close all
        end
    end
    
end

%% computing quality of rampiness - pearson/spearman and significance of ramp

options.ind_train = screeningData.imageIDs; % use all objects to calculate STA
screeningData.sigRampSTA = ones(length(strctCells), length(firstWindow:stepSize:lastWindow));
screeningData.corrRatioRamp = ones(length(strctCells), length(firstWindow:stepSize:lastWindow))*-1;
% significance of ramps 
tic
alpha = 0.01;
for cellIndex = l(strctCells)
    if ~isempty(screeningData.responses{cellIndex, 2}) % responsive cells only
        
        for bin = 1:length(firstWindow:stepSize:lastWindow)
            
            resp = screeningData.rollingWindowSTA{cellIndex, 1}(:, bin);
            [p, corrRatio] = linearity_measure_STA(resp, params, options);
            screeningData.sigRampSTA(cellIndex, bin) = p;
            screeningData.corrRatioRamp(cellIndex, bin) = corrRatio;
            
        end
        
        windowRange = firstWindow+(screeningData.timelimits(1)*1e3):stepSize:lastWindow+(screeningData.timelimits(1)*1e3);
        yc_ev = find(windowRange == screeningData.responses{cellIndex, 2});
        isSig = find(screeningData.sigRampSTA(cellIndex, :) < alpha);
        if ~isempty(isSig)
            % plot relationship
            f = figure; clf;
            hold on
            plot(firstWindow+(screeningData.timelimits(1)*1e3):stepSize:lastWindow+(screeningData.timelimits(1)*1e3), screeningData.sigRampSTA(cellIndex, :));
            set(gca,'FontWeight', 'bold')
            yline(alpha, '--',  'FontWeight', 'bold');
            title({['Sig and Correlation Ratio vs Rolling Window'], [strctCells(cellIndex).brainArea '\_' num2str(strctCells(cellIndex).ChannelNumber) '\_' num2str(strctCells(cellIndex).Name)]});
            xlabel('WindowBegin post stimON', 'FontWeight', 'bold');
            ylabel('p-value of ramp', 'FontWeight', 'bold');
            
            yyaxis right
            plot(firstWindow+(screeningData.timelimits(1)*1e3):stepSize:lastWindow+(screeningData.timelimits(1)*1e3), screeningData.corrRatioRamp(cellIndex, :));%, 'Color', [0.4660 0.6740 0.1880]);
            yline(1, '--',  'FontWeight', 'bold');
            ylabel('Pearson\Spearman', 'FontWeight', 'bold');
            
            yyaxis left
            if ~isempty(yc_ev)
                plot(screeningData.responses{cellIndex, 2}, screeningData.sigRampSTA(cellIndex, yc_ev), '*','MarkerSize',10, 'MarkerEdgeColor', [0.4940 0.1840 0.5560]);
            end
            
            yyaxis right
            if ~isempty(yc_ev)
                plot(screeningData.responses{cellIndex, 2}, screeningData.corrRatioRamp(cellIndex, yc_ev), '*','MarkerSize',10);%, 'MarkerEdgeColor', [0.4940 0.1840 0.5560]);
            end
            
            print(f, [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_Linearity'], '-dpng', '-r0');
            close all
        end
    end
end
toc

%% Plot cosine similarity between the axis I ended up choosing and all the time windows




















