
setDiskPaths

load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create dataParams = 500x50

% task = 'Object_Screening';
% binnedAverage = 1;

task = 'Recall_Task';
binnedAverage = 0;

separateSides = 0;

%% choose task to analyze and load in data accordingly

taskPath = [diskPath filesep task];

if strcmp(task, 'Object_Screening')
    
    load([taskPath filesep 'MergedITCells_500stim_Scrn_SigRamp']); % loads cells with sig ramps
    taskPath = [taskPath filesep 'forPaper'];
    if ~exist(taskPath)
        mkdir([taskPath]);
    end

elseif strcmp(task, 'Recall_Task')  
    
%     load([diskPath filesep 'ObjectSpace' filesep 'parameters_500_objects.mat']); % will create dataParams = 500x50
    load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500stim_Im_SigRamp']);
    load([diskPath filesep 'Recall_Task' filesep 'AllITResponses_500stim_Im_SigRamp']);
%     load([diskPath filesep 'Recall_Task' filesep 'AllCells_500stim_Im']);
%     load([diskPath filesep 'Recall_Task' filesep 'AllResponses_500stim_Im']);
    taskPath = [taskPath filesep 'forPaper'];
    if ~exist(taskPath)
        mkdir([taskPath]);
    end
end
strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';


% cleave out garbage neurons - improves things
% only 3 of them are in sigramp set
% badCells = [1085, 1759, 713, 5477, 2099, 750, 773, 819]'; 
% allCells = cell2mat(strctCELL(:, 1));
% goodCells = ~ismember(allCells, badCells);
% 
% strctCells = strctCells(goodCells);
% strctCELL = strctCELL(goodCells, :);
% strctResp = strctResp(goodCells);
% psths = psths(goodCells, :);
% responses = responses(goodCells, :);


%% separating hemispheres
if separateSides
    left = 0;
    right = 1;
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    if left
        ITSide = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));
        
    elseif right
        ITSide = cellfun(@(x) strcmp(x, 'RFFA'), strctCELL(:, 4));
        
    end
    
    strctCells = strctCells(ITSide);
    psths = psths(ITSide, :);
    responses = responses(ITSide, :);
    if strcmp(task, 'Recall_Task')
        strctResp = strctResp(ITSide);
    end
end
%% compute binned averages
numOfAxes = 2;

distMat_pref = zeros(length(strctCells), 8);
distMat_ortho = zeros(length(strctCells), 8);

if binnedAverage
    
    image_IDs = [1:500]';
    options.ind_train = image_IDs; % regardless of task (screening/recall) you need all 500 imgs to compute axis
    options.screenType = 'Object';
    options.axis = 'STA';
    if strcmp(task, 'Object_Screening')
        options.nbins = 20;
    elseif strcmp(task, 'Recall_Task')
        options.nbins = 5;
    end
    % firing_rate = cell(length(strctCells), 1);
    f_r = zeros(length(strctCells), options.nbins);
    f_r_o = zeros(length(strctCells), options.nbins);
    max_f_r = zeros(length(strctCells), 1);
    
    for cellIndex = 1:length(strctCells)
        if strcmp(task, 'Recall_Task')
            options.recalled_stim = strctCells(cellIndex).recalledStim;
        end
        nonlin = Utilities.ObjectSpace.binned_average_perCell_STA(responses{cellIndex, 1}, params, options);
        %     firing_rate{cellIndex, 1} = nonlin.y./max(abs(nonlin.y));
        max_f_r(cellIndex) = max(abs(nonlin.y));
        f_r(cellIndex, :) = (nonlin.y./max_f_r(cellIndex))';
        
        % getting rid of nans - does this make sense??
        if ~isempty(find(isnan(f_r(cellIndex, :))))
            ids = find(isnan(f_r(cellIndex, :)));
            f_r(cellIndex, ids) = 0;
            %         ids = setdiff(ids, 1);
            %         if ~isempty(ids)
            %             lid = min(ids);
            %             rid = max(ids);
            %             if lid ~= 1 && rid ~= options.nbins
            %                 interp = linspace(f_r(cellIndex, lid-1), f_r(cellIndex, rid+1), length(ids)+2);
            %                 f_r(cellIndex, ids) = interp(2:end-1);
            %                 %         else
            %             end
            %         end
        end
        
    end
    
    
    
    options.axis = 'ortho';
    
    for cellIndex = 1:length(strctCells)
        if strcmp(task, 'Recall_Task')
            options.recalled_stim = strctCells(cellIndex).recalledStim;
        end
        nonlin = Utilities.ObjectSpace.binned_average_perCell_STA(responses{cellIndex, 1}, params, options);
        
        %     firing_rate{cellIndex, 2} = nonlin.y./max(abs(nonlin.y));
        f_r_o(cellIndex, :) = (nonlin.y./max_f_r(cellIndex))';
        
        if ~isempty(find(isnan(f_r_o(cellIndex, :))))
            ids = find(isnan(f_r_o(cellIndex, :)));
            f_r_o(cellIndex, ids) = 0;
        end
    end
elseif ~binnedAverage && strcmp(task, 'Recall_Task')
    % for recall only - not binned average, line up responses along axis

    for ax = 1:numOfAxes
        if ax == 1
            axis = 'STA';
        elseif ax == 2           
            axis = 'ortho';
        end
        
        if strcmp(axis, 'STA')
            f_r = zeros(length(strctCells), 8);
            max_f_r = zeros(length(strctCells), 1);
            f_r_hist = zeros(length(strctCells), 8); % for histogram of slopes
            
        elseif strcmp(axis, 'ortho')
            f_r_o = zeros(length(strctCells), 8);
            f_r_o_hist = zeros(length(strctCells), 8);

        end
        
        for cellIndex = 1:length(strctCells)
            if strcmp(axis, 'STA')
                crds = strctCells(cellIndex).Im_xvals;
                
                [~, idx] = sort(crds);
                resp = strctResp(cellIndex).CRResp';
                resp = resp(idx);
                raw_resp = resp;
                resp = resp - mean(resp);
                max_f_r(cellIndex) = max(resp);
                
                
                if length(strctCells(cellIndex).recalledStim) == 6 % first 9 cells only imagined 6 stim
                    f_r(cellIndex, 3:end) = resp./max_f_r(cellIndex);
                    f_r_hist(cellIndex, 3:end) = resp./std(resp);  % fpr histogram don't want to divide by max                    
                    distMat_pref(cellIndex, 3:end) = crds(idx)./max(crds);
                    
                else
                    f_r(cellIndex, :) = resp./max_f_r(cellIndex);
                    f_r_hist(cellIndex, :) = resp./std(resp);  
                    distMat_pref(cellIndex, :) = crds(idx)./max(crds);
                end
                
            elseif strcmp(axis, 'ortho')
                assert(length(unique(max_f_r)) > 1, 'Need to compute max for normalization');
                crds = strctCells(cellIndex).Im_yvals';
                [~, idx] = sort(crds);
                resp = strctResp(cellIndex).CRResp';
                resp = resp(idx);
                raw_resp = resp;
                resp = resp - mean(resp);
                
                if length(strctCells(cellIndex).recalledStim) == 6
                    f_r_o(cellIndex, 3:end) = resp./max_f_r(cellIndex); % using common scale
                    f_r_o_hist(cellIndex, 3:end) = resp./std(resp); % not the standardization here
                    distMat_ortho(cellIndex, 3:end) = crds(idx)./max(crds);
                else
                    f_r_o(cellIndex, :) = resp./max_f_r(cellIndex);
                    f_r_o_hist(cellIndex, :) = resp./std(resp);
                    distMat_ortho(cellIndex, :) = crds(idx)./max(crds);
                end
                
            end
            
        end
    end
end

%%
% smoothing for plotting heatplot
f_r_s = f_r;
f_r_os = f_r_o;

if exist('f_r_hist', 'var')
    f_r_s_hist = f_r_hist;
    f_r_os_hist = f_r_o_hist;
end
binsize = 2;
for cellIndex = 1:length(strctCells)
    f_r_s(cellIndex, :) = Utilities.Smoothing.fastsmooth(f_r(cellIndex, :), binsize, 1, 1);
    f_r_os(cellIndex, :) = Utilities.Smoothing.fastsmooth(f_r_o(cellIndex, :), binsize, 1, 1);
    
    if exist('f_r_hist', 'var')
        
        f_r_s_hist(cellIndex, :) = Utilities.Smoothing.fastsmooth(f_r_hist(cellIndex, :), binsize, 1, 1);
        f_r_os_hist(cellIndex, :) = Utilities.Smoothing.fastsmooth(f_r_o_hist(cellIndex, :), binsize, 1, 1);
    end
end

meanfro = mean(f_r_os, 1);
meanfr = mean(f_r_s, 1);

% meanfro = mean(f_r_o, 1);
% meanfr = mean(f_r, 1);


%% make heat plot - make this a function 
useSmoothed = 1;
avgCells = 1;
% ax2plot = 'sta';
ax2plot = 'ortho';


f = figure; 
% hold on
if avgCells
    if useSmoothed
        meanfro = mean(f_r_os, 1);
        meanfr = mean(f_r_s, 1);
    else
        meanfro = mean(f_r_o, 1);
        meanfr = mean(f_r, 1);
    end
    mat = [meanfr; meanfro];

else
    if useSmoothed
        if strcmp(ax2plot, 'sta')
            mat = f_r_s;
        elseif strcmp(ax2plot, 'ortho')
            mat = f_r_os;
        end
        sortBy = f_r_s;
    else
        if strcmp(ax2plot, 'sta')
            mat = f_r;
        elseif strcmp(ax2plot, 'ortho')
            mat = f_r_o;
        end
        sortBy = f_r;
    end
    
    [~, b] = sortrows(sortBy(:, end), 'descend');
    mat = mat(b, :);

end

% set filename ------------------------------------------------------------
if useSmoothed
    titleSmth = ['_smoothed' num2str(binsize)];
else
    titleSmth = '';
end
if strcmp(task, 'Recall_Task')
    if avgCells
        filename = [taskPath filesep 'Im_RampSummary_AcrossCells'  titleSmth];
    else
        filename = [taskPath filesep 'Im_RampSummary_' ax2plot  titleSmth];        
    end
elseif strcmp(task, 'Object_Screening')
    if avgCells
        filename = [taskPath filesep 'Scrn_RampSummary_AcrossCells' titleSmth];
    else
        filename = [taskPath filesep 'Scrn_RampSummary_' ax2plot '_' num2str(options.nbins) titleSmth];
    end
end

% Plot map and set colors ----------------------------------------
% imagesc(mat, [0 max(max(mat))]); 
% imagesc(mat, [0 1]); 
% f = imagesc(f_r, [-1 1]);
imagesc(mat, [-max(max(mat)) max(max(mat))]); 

orangemap = esa(300);
[WhiteColor, Whitepos] = max(sum(orangemap,2));
orangemap = orangemap([1:Whitepos round(linspace(Whitepos+1,size(orangemap,1)-2,Whitepos-1))],:);
colormap(orangemap)

% set up colorbar -------------------------------------------------
cb = colorbar;
cb.Ticks = ([-2 2]);
cb.FontWeight = 'bold';
cb.Label.String = 'Norm Resp';
cb.Position = [0.920833333333333,0.719047619047619,0.025595238095239,0.204761904761913];

% set title and axis labels correctly -----------------------------------
if strcmp(task, 'Recall_Task')
    taskTitle = 'Imagination';
elseif strcmp(task, 'Object_Screening')
    taskTitle = 'Screening';
end

if avgCells 
    tasksubTitle = 'average across cells';
    yticklabels({'', 'STA', '', 'ortho', ''});
    ylabel('Average Across Condition', 'FontWeight', 'bold');

else
    ylabel('Cell number', 'FontWeight', 'bold');
    if strcmp(ax2plot, 'sta')
        tasksubTitle = 'distance along preferred axis';
    elseif strcmp(ax2plot, 'ortho')
        tasksubTitle = 'distance along principal orthogonal axis';
    end
end
title([taskTitle '-' tasksubTitle], 'FontSize', 14, 'FontWeight', 'bold');
if strcmp(task, 'Recall_Task')
%     if avgCells
%         xlabel('Normalized Distance', 'FontWeight', 'bold');
%         set(gca, 'xtick', [0.5 size(f_r, 2)+0.5]);
%         set(gca, 'xticklabel', [-1 1], 'FontWeight', 'bold');
%     else
        xlabel('Stimuli (sorted by distance)', 'FontWeight', 'bold');
%     end
elseif strcmp(task, 'Object_Screening')
    xlabel('Normalized Distance', 'FontWeight', 'bold');
    set(gca, 'xtick', [0.5 size(f_r, 2)+0.5]);
    set(gca, 'xticklabel', [-1 1], 'FontWeight', 'bold');
end
set(gca, 'FontWeight', 'bold');
if exist('left', 'var') && left
    filename = [filename '_LeftIT'];

elseif exist('right', 'var') && right
    filename = [filename '_RightIT'];

end
filename = [filename '_forPaper'];

% print(f, filename, '-dpng', '-r0');
% close all

%% histogram of slopes 

% take responses to im stimuli sorted by pref and ortho ax distance computed above (f_r and f_r_o)
% plot them as scatter per cell
% fit a line through them and record slope
% plot histogram

cols = Utilities.distinguishable_colors(length(strctCells));
dot_size = 20;
slopes_pref = zeros(length(strctCells), 1);
slopes_ortho = zeros(length(strctCells), 1);

cc = zeros(length(strctCells), 2);

R_Squared = zeros(length(strctCells), 2);
for ax = 1:numOfAxes
    
    if ax == 1
        matToUse = f_r_hist;
        distMat = distMat_pref;
%         matToUse = f_r;
    elseif ax == 2
        matToUse = f_r_o_hist; 
        distMat = distMat_ortho;
%         matToUse = f_r_o;       
    end
    
%     f = figure;
%     hold on
    for cellIndex = 1:length(strctCells)
        
        c_x = distMat(cellIndex, end-(length(strctCells(cellIndex).recalledStim)-1):end);
        x = 1:length(strctCells(cellIndex).recalledStim);
        y = matToUse(cellIndex, end-(length(strctCells(cellIndex).recalledStim)-1):end);
        % doing this strange end-n bullshit because not all rows are of the same length
%         scatter(x, y, dot_size, cols(cellIndex,:), 'filled')
        
        % best fit line
        % Get coefficients of a line fit through the data.
        [coefficients, stats] = polyfit(c_x, y, 1);
        
        cc(cellIndex, ax) = corr(c_x', y');
        
        if ax == 1            
            slopes_pref(cellIndex) = coefficients(1);
        elseif ax == 2
            slopes_ortho(cellIndex) = coefficients(1);

        end
        
        % Create a new x axis with exactly 1000 points (or whatever you want).
        xFit = linspace(0, 8, 1000);
        % Get the estimated yFit value for each of those 1000 new x locations.
        yFit = polyval(coefficients , xFit);
        
        R_Squared(cellIndex, ax) = 1 - (stats.normr/norm(y - mean(y)))^2;

        %         plot(xFit, yFit, 'Color', cols(cellIndex, :), 'LineWidth', 1.5); % Plot fitted line.

        
        
    end
    
end



%% trying to make mirrored histogram
% histogram
% -0.3:0.04:0.4 binedges best vis so far

% test if distributions are significantly different 
[h, p] = kstest2(slopes_pref, slopes_ortho);
 mirrored = 0;

f = figure; 
hold on; 

if ~mirrored
%     histogram(slopes_ortho, 'BinEdges', -0.3:0.04:0.4)
%     histogram(slopes_pref, 'BinEdges', -0.3:0.04:0.4)
    histogram(slopes_ortho)
    histogram(slopes_pref)
    x_lim = xlim;
    y_lim = ylim;
    filename = [taskPath filesep 'SlopesComaprison_BestFitLines_STAvsOrtho'];

else
    h1 = axes;
    histogram(h1, slopes_ortho, 'BinEdges', -0.3:0.04:0.4, 'FaceColor', [0 0.4470 0.7410])
    set(h1, 'YDir', 'reverse')
    set(h1, 'YAxisLocation', 'Right')
    set(h1, 'Xtick', []);
    set(h1, 'Ytick', []);
    
    ylim([-10 10])
    
    h2 = axes;
    histogram(h2, slopes_pref, 'BinEdges', -0.3:0.04:0.4, 'FaceColor', [0.4940 0.1840 0.5560])
    % set(h2, 'XLim', get(h1,'Xlim'))
    set(h2, 'Color', 'None')
    % set(h2, 'Xtick', []);
    ylim([-10 10])
    x_lim = xlim;
    
    filename = [taskPath filesep 'SlopesComaprison_BestFitLines_STAvsOrtho_mirrored'];

end
text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');
lgnd = legend({'Ortho axis', 'Preferred axis'});
xlabel('Slope of best fit line');
ylabel('No of neurons');
title({['Slopes of best fit lines - STA and Orthogonal axes'], 'Standardized FR'})
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% print(f, filename, '-dpng', '-r0')
% close all




%% correlation histogram
f = figure; 
hold on; 
histogram(cc(:, 2), 'BinEdges', -1:0.1:1) % ortho
histogram(cc(:, 1), 'BinEdges', -1:0.1:1) % pref

% test if distributions are significantly different 
[h, p] = kstest2(cc(:, 1), cc(:, 2));


x_lim = xlim;
% ylim([0 15])
y_lim = ylim;
text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');
lgnd = legend({'Ortho axis', 'Preferred axis'});
xlabel('Correlation value');
ylabel('No of neurons');
title({'Correlation of projection value vs firing rate', 'STA and Orthogonal axes'})
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
filename = [taskPath filesep 'CorrProjNormValvsFR_STAvsOrtho'];
% print(f, filename, '-dpng', '-r0')
% close all


%% cdfs - Slopes 

e_cdf = 0;
f = figure; 
hold on; 

if e_cdf
    % % ecdf with bounds - try to make this work
    ecdf(slopes_ortho, 'Bounds', 'on'); % ortho
    ecdf(slopes_pref, 'Bounds', 'on');% pref
    grid on
    % plot(f1, x1, 'LineWidth', 2);
    % plot(f2, x2, 'LineWidth', 2);
    
else
    cd1 = cdfplot(slopes_ortho);
    cd2 = cdfplot(slopes_pref);
    cd1.LineWidth = 2;
    cd2.LineWidth = 2;
end
% test if distributions are significantly different 
[h, p] = kstest2(slopes_ortho, slopes_pref);


x_lim = xlim;
% ylim([0 15])
y_lim = ylim;
text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');

if e_cdf
    lgnd = legend({'Ortho axis','', '', 'Preferred axis','', ''});
    filename = [taskPath filesep 'ECDFBestFitSlope_STAvsOrtho'];

else
    lgnd = legend({'Ortho axis','Preferred axis'});
    filename = [taskPath filesep 'CDFBestFitSlope_STAvsOrtho'];

end
lgnd.Position = [0.200892857142857,0.684523809523811,0.224107142857143,0.08452380952381];
xlabel('x = Slope of best fit line');
% ylabel('No of neurons');
title({'Correlation of projection value vs firing rate', 'STA and Orthogonal axes'})
set(gca, 'FontSize', 14, 'FontWeight', 'bold');


% print(f, filename, '-dpng', '-r0')



%% cdfs - correlation

e_cdf = 0;
f = figure; 
hold on; 

if e_cdf
    % % ecdf with bounds - try to make this work
    e1 = ecdf(cc(:, 2), 'Bounds', 'on'); % ortho
    e2 = ecdf(cc(:, 1), 'Bounds', 'on');% pref
    grid on
%     plot(e1, x1, 'LineWidth', 2);
%     plot(e2, x2, 'LineWidth', 2);
    
else
    cd1 = cdfplot(cc(:, 2));
    cd2 = cdfplot(cc(:, 1));
    cd1.LineWidth = 2;
    cd2.LineWidth = 2;

    [c_p, x_p, ~, ~, ~] = cdfcalc(cc(:, 1));
    [c_o, x_o, ~, ~, ~] = cdfcalc(cc(:, 2));
end
% test if distributions are significantly different 
[h, p] = kstest2(cc(:, 1), cc(:, 2));


x_lim = xlim;
% ylim([0 15])
y_lim = ylim;
text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');

if e_cdf
    lgnd = legend({'Ortho axis','', '', 'Preferred axis','', ''});
    filename = [taskPath filesep 'ECDFProjvsFR_STAvsOrtho'];

else
    lgnd = legend({'Ortho axis','Preferred axis'});
    filename = [taskPath filesep 'CDFProjvsFR_STAvsOrtho'];

end
lgnd.Position = [0.200892857142857,0.684523809523811,0.224107142857143,0.08452380952381];
xlabel('x = Correlation value');
% ylabel('No of neurons');
title({'Correlation of projection value vs firing rate', 'STA and Orthogonal axes'})
set(gca, 'FontSize', 14, 'FontWeight', 'bold');


% print(f, filename, '-dpng', '-r0')






