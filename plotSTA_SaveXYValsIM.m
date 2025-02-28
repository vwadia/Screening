% Script to load in all cells and plot STA plots for all of them
% If recall sessions - mark recalled stim
% 
% Note: Normally you do this 1 session at a time in screeningScript - but after large IT Resort I need to do this


%% set paths and load data

setDiskPaths

% taskPath = 'Object_Screening';
% markPoints = false;
% matchColsToRecall = false;

taskPath = 'Recall_Task';
markPoints = false;
matchColsToRecall = false; % this doesn't work quite yet
% [sessID, ~, recStim] = Utilities.sessionListAllTasks(taskPath, false, true);


options.noTicks = false;

if strcmp(taskPath, 'Object_Screening')
    load([diskPath filesep taskPath filesep 'MergedITCells_500Stim_Scrn_SigRamp'])
elseif strcmp(taskPath, 'Recall_Task')
    load([diskPath filesep taskPath filesep 'AllITCells_500Stim_Im_SigRamp'])
end

pathOut = [diskPath filesep taskPath filesep 'forPaper' filesep 'STA_and_projections_varWindow_sigCells'];


if isfield(options, 'noTicks') && options.noTicks
    pathOut = [pathOut '_noAxTicks'];
    
    if markPoints
        pathOut = [pathOut '_markedPoints'];
    end
end




%%
% layermat = 'fc6'; stimDir = '500Stimuli';
% load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_AlexnetPYTHON_MatlabMean_' layermat '_' stimDir '.mat']]); params = feat;

load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create params = 500x50

if size(params, 2) == 4096
    pathOut = [pathOut '_4096D'];
end

imageIDs = [1:500]';
options.screenType = 'Object';

% options.recalledCols = [1 0.5 0;...% orange
%     1 0.25 0;...% orange
%     0.8 1 0;...% yellow
%     0 1 0;...% green
%     0 1 1;...% turquoise
%     1 0 1;...% pink
%     0.75 0.75 0.75;...% grey
%     0 0 0];% white 

%%

if matchColsToRecall
    pathOut = [pathOut '_ColsMatchedtoRecall'];
end

if ~exist(pathOut)
    mkdir([pathOut]);
end
if ~exist('sessID', 'var')
    sessID = 1;
end
    
for ss = 1:length(sessID)
    
    
    
    
for cellIndex = 1:length(strctCells)

    
    options.ind_train = imageIDs; % use all objects to calculate STA

    if strcmp(taskPath, 'Recall_Task')
        if markPoints
            if strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen') || strcmp(strctCells(cellIndex).SessionID, 'P76CSFast_2')
                options.marked_positions = [1 456 457 221 376 114]; % P76 Recall 1
                options.recalled_stim = [12 19 25 123 270 487]; % P76 Recall 1
            elseif strcmp(strctCells(cellIndex).SessionID, 'P76CS_RecScreen3')
                options.marked_positions = [16 34 144 358 382 450]; % P76 Recall 2
                options.recalled_stim = [54 129 130 186 270 449]; % P76 Recall 2
            elseif  strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen_3') || strcmp(strctCells(cellIndex).SessionID, 'P76CS_RecScreen_3')
                options.marked_positions = [175 107 345 340 476 459 499 496]; % P76 Recall 3
                options.recalled_stim = [18 44 45 81 135 181 230 344];
            elseif  strcmp(strctCells(cellIndex).SessionID, 'P79CS_1') || strcmp(strctCells(cellIndex).SessionID, 'P79CS_ReScreen_1')
                options.marked_positions = [3 104 167 196 422 453 473 491];
                options.recalled_stim = [9 157 167 200 201 291 422 498];
            elseif  strcmp(strctCells(cellIndex).SessionID, 'P79CS_3') || strcmp(strctCells(cellIndex).SessionID, 'P79CS_ReScreen_3')
                options.marked_positions = [476 497 453 495 279 309 13 167];
                options.recalled_stim = [9 12 117 292 360 368 421 492];
            elseif  strcmp(strctCells(cellIndex).SessionID, 'P79CS_4') || strcmp(strctCells(cellIndex).SessionID, 'P79CS_ReScreen_4')
                options.marked_positions = [49 60 70 173 422 480 499 500];
                options.recalled_stim = [77 112 160 232 278 345 387 440];
            elseif  strcmp(strctCells(cellIndex).SessionID, 'P80CS_RecScreen_1') || strcmp(strctCells(cellIndex).SessionID, 'P80CS_ReScreenRecall')
                options.marked_positions = [ 5   185   388   431   462   468   486   498];
                options.recalled_stim = [17    61    76   114   157   161   177   480 ];
            elseif strcmp(strctCells(cellIndex).SessionID, 'P80CS_RecScreen_2') || strcmp(strctCells(cellIndex).SessionID, 'P80CS_ReScreecRecall_2')
                options.marked_positions = [134 257 454 466 488 492 498 499];
                options.recalled_stim = [55 88 148 251 256 274 285 365];
            elseif strcmp(strctCells(cellIndex).SessionID, 'P84CS_RecScreen_1') || strcmp(strctCells(cellIndex).SessionID, 'P84CS_ReScreenRecall_1')
                options.marked_positions = [144   175   260   281   325   387   466   478];
                options.recalled_stim = [68 121 243 261 281 308 415 434];
            elseif strcmp(strctCells(cellIndex).SessionID, 'P84CS_RecScreen_2') || strcmp(strctCells(cellIndex).SessionID, 'P84CS_ReScreenRecall_2')
                options.marked_positions = [18   191   297   329   379   466   471   483];
                options.recalled_stim = [52 141 195 223 238 351 399 482];
            elseif strcmp(strctCells(cellIndex).SessionID, 'P85CS_RecScreen_1') || strcmp(strctCells(cellIndex).SessionID, 'P85CS_ReScreenRecall')
                options.marked_positions = [ 4   160   306   329   387   397   497   500];
                options.recalled_stim = [7 11 99 153 201 251 355 486];
            end
        end
        
        
        if ~isfield(options, 'recalled_stim')
            
            matchColsToRecall = false;
        end
        if matchColsToRecall
            options.recalledCols = Utilities.distinguishable_colors(length(options.recalled_stim));
            
        elseif ~matchColsToRecall && isfield(options, 'recalled_stim')
            
            if length(options.recalled_stim) == 8
                options.recalledCols = [1 0.5 0;...% orange
                    1 0.25 0;...% orange
                    0.8 1 0;...% yellow
                    0 1 0;...% green
                    0 1 1;...% turquoise
                    1 0 1;...% pink
                    0.75 0.75 0.75;...% grey
                    0 0 0];% white
                
            elseif length(options.recalled_stim) == 6
                options.recalledCols = [1 0.5 0;...% orange
                    0.8 1 0;...% yellow
                    0 1 0;...% green
                    0 1 1;...% turquoise
                    1 0 1;...% pink
                    0.75 0.75 0.75];% grey
                
            end
        end
        options.matchColsToRecall = matchColsToRecall;
        
        %% HAVE TO FIX STRCTRESP BEFORE I CAN WRITE THIS IN THIS WAY
        if matchColsToRecall
            % have to reset colors to shuffling is correct
            options.recalledCols = Utilities.distinguishable_colors(length(options.recalled_stim));
            options.recalledCols = circshift(options.recalledCols, -1); % to keep colors consistent with my stupid oroginal way of plotting

            load([basePath filesep 'ITResponses.mat'])
            options.task = 'Recall_Task';
            if isequal(screeningData.responses{cellIndex, 3}, strctResp(resp_ctr).Name)
                options.recalledStim = options.recalled_stim;
                options.ScrnResp = strctResp(resp_ctr).ScrnResp;
                options.CRResp = strctResp(resp_ctr).CRResp;
                
                options.cellName = strctResp(resp_ctr).Name;
            else
                keyboard
            end
            % important bit - figure out order along axis
            axToUse = 'sta';
            orderAlongAx = Utilities.ObjectSpace.returnOrderAlongAx(screeningData.responses{cellIndex, 1}, params, options, axToUse);
           
            options.recalledCols = options.recalledCols(orderAlongAx, :);
            options.axToUse = axToUse;
        end
        
        
     end
    
    
    
    [hfig, p, options] = Utilities.ObjectSpace.STA_figure_original(responses{cellIndex, 1}, params, options); % pass score to this instead of projectedResponses
    if isfield(options, 'xvals') && isfield(options, 'yvals')
        strctCells(cellIndex).Im_xvals = options.xvals;
        strctCells(cellIndex).Im_yvals = options.yvals;
    end
    
    if isfield(options, 'recalled_stim')
        strctCells(cellIndex).recalledStim = options.recalled_stim;
    end
    
    
%     sgtitle({['Cell number ' num2str(strctCells(cellIndex).Name)] 'STA and projections ', [strctCells(cellIndex).brainArea ' - ' strctCells(cellIndex).SessionID]});
%     sgtitle({['Cell number ' num2str(strctCells(cellIndex).Name)] 'STA and projections ', strctCells(cellIndex).brainArea});
    if isfield(options, 'noTicks') && options.noTicks
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.2);
        %         set(findobj(gcf,'type','axes'),'FontName','Arial','FontWeight','Bold', 'LineWidth', 1.2);
    else
%         sgtitle({['Cell number ' num2str(strctCells(cellIndex).Name)] 'STA and projections ', strctCells(cellIndex).brainArea});
            set(findobj(gcf,'type','axes'),'FontName','Arial','FontWeight','Bold', 'LineWidth', 1.2);

    end

%     set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    
    if isfield(options, 'encoded_stim') || isfield(options, 'recalled_stim')
        filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_orderedStim_' strctCells(cellIndex).SessionID];
    else
        filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) '_' strctCells(cellIndex).SessionID];

    end
    
    print(hfig, filename, '-dpng', '-r0')

    
    close all
    
    
    

end
end
% save([diskPath filesep taskPath filesep 'AllITCells_500Stim_Im_SigRamp'], 'strctCells', 'responses', 'psths', '-v7.3')

%%
% p = [];
% for cellIndex = l(strctCells)
%     
%     
%     [p(cellIndex), ~] = Utilities.ObjectSpace.linearity_measure_STA(responses{cellIndex, 1}, params, options);
%     
%     
%     
%     
% end



