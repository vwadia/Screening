% Script to load in synthetic psths and responses
% and see if the FRs increase to the selected images in gradual steps
% 
% Note:
% The cells must already be merged 
% The list of chosen morning cells must be inputted manually
% 
% Load in synth responses (hopefully these are sorted correctly) then  
% plot a scatter plot of them

setDiskPaths

patID = 'P81CS';
patID = 'P82CS';


if strcmp(patID, 'P81CS')
    a_basePath = [diskPath filesep 'Object_Screening' filesep patID filesep 'ClosedLoopReScreen_Session_1_20221030'];
    m_basePath = [diskPath filesep 'Object_Screening' filesep patID filesep 'ClosedLoopScreening_Session_1_20221030'];

elseif strcmp(patID, 'P82CS')
    a_basePath = [diskPath filesep 'Object_Screening' filesep patID filesep 'ClosedLoopReScreen_Session_1_20230115'];
    m_basePath = [diskPath filesep 'Object_Screening' filesep patID filesep 'ClosedLoopScreening_Session_1_20230115'];

end

load([a_basePath filesep 'SynthPsthandResponses.mat'])
load([a_basePath filesep 'strctCells.mat'])

%% computing correlation between morn and aft axes


morn = load([m_basePath filesep 'PsthandResponses']);
aft = load([a_basePath filesep 'PsthandResponses']);

layermat = 'fc6';

stimDir = '500Stimuli';
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);

ax_morn = [];
ax_aft = [];

if strcmp(patID, 'P81CS')
    corrThresh = 0.9;
    morn_cells = [1100; 1546; 1877; 1897; 2156];
    aft_matches = [1525; 1943; 1928; 2300; 1935];
    aft_idx = [1; 2; 3; 4; 5];  
    legNames = {'1525', '1943', '1928', '2300', '1935'};

elseif strcmp(patID, 'P82CS')
    corrThresh = 0.9;
    morn_cells = [904];
    aft_matches = [2360]; 
    aft_idx = [1];
    legNames = {'2360'};

end


for idx = 1:length(aft_matches)
    
    cellId = cellfun(@(x) isequal(x, aft_matches(idx)), aft.responses(:, 3), 'UniformOutput', false);
    
    cellIndex = find(cell2mat(cellId) == 1);
    
    [ax_aft(idx, :), ~] = Utilities.ObjectSpace.analysis_STA(aft.responses{cellIndex, 1}, params, 'sta');
    
end



% compute axes for all morning cells that are matched to afternoon
for idx = 1:length(aft_matches)
    
    cellId = cellfun(@(x) isequal(x, morn_cells(aft_idx(idx))), morn.responses(:, 3), 'UniformOutput', false);
    
    cellIndex = find(cell2mat(cellId) == 1);
    
    [ax_morn(idx, :), ~] = Utilities.ObjectSpace.analysis_STA(morn.responses{cellIndex, 1}, params, 'sta');

    
end

axcorr = diag(corr(ax_morn', ax_aft'));

cells2use = axcorr > corrThresh;
aft_matches = aft_matches(cells2use);
aft_idx = find(cells2use == 1);

if strcmp(patID, 'P81CS')
    num_steps = 20; % 10 in each direction - inputed manually
    num_reps = 4;
    num_trials = num_steps*4;
    scale = 1;
elseif strcmp(patID, 'P82CS')
    num_steps = 10; % inputed manually
    num_reps = 4;
    num_trials = num_steps*4;
    scale = 2;
end

% cols  = Utilities.distinguishable_colors(length(aft_matches));
cols = [0 0.4470 0.7410;...	
0.8500 0.3250 0.0980;...	
0.9290 0.6940 0.1250;...
0.4940 0.1840 0.5560;...	
0.4660 0.6740 0.1880;...	
0.3010 0.7450 0.9330;...
0.6350 0.0780 0.1840];	

%% plotting all cells together

f1 = figure; 
hold on
title({'Response to synthetic images', ['Mean axis corr ' num2str(mean(axcorr(aft_idx)))]})

for cellIndex = 1:length(aft_matches)


    cellId = cellfun(@(x) isequal(x, aft_matches(cellIndex)), synthResponses(:, 3), 'UniformOutput', false);
    
    id = find(cell2mat(cellId) == 1);
    
    cellName = cell2mat(synthResponses(cell2mat(cellId), 3))
    
    
    windowBegin = synthResponses{id, 2} + 170;
    windowEnd = windowBegin+267;
    
    cellPsth = synthPsth{id, 1}(:, windowBegin:windowEnd);
    
    
    
    Resp = mean(cellPsth, 2)*1e3;
    cellResp = Resp((aft_idx(cellIndex)-1)*num_trials+1:aft_idx(cellIndex)*num_trials, 1);
    
    cellResp = reshape(cellResp, [num_reps num_steps]);
    
    
    err = std(cellResp)./size(cellResp, 1);
    e = errorbar(1:num_steps, mean(cellResp, 1), err, 'LineWidth', 2);
    e.Color = cols(cellIndex, :);
    e.Marker = 'o';
    e.MarkerFaceColor = cols(cellIndex, :);
    e.MarkerEdgeColor = cols(cellIndex, :);
    
%     cellName = cell2mat(synthResponses(cell2mat(cellId), 3));
%     
%     cellResp = synthResponses{cell2mat(cellId), 1};
%     
%     p = plot(1:num_steps, cellResp((aft_idx(cellIndex)-1)*num_steps+1:aft_idx(cellIndex)*num_steps, 1), '-o')%, 'Color', cols(cellIndex, :))
%     p.LineWidth = 2;
    

    
end
xtiklabs = {setdiff(-num_steps/2:scale:num_steps/2, 0)};
xticks([1:20])
xticklabels(xtiklabs)
xlabel('Steps along preferred axis')
ylabel('Spikes/s')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
lgnd = legend(legNames(aft_idx), 'FontWeight', 'bold');

if ~exist('cells2use', 'var')
    lgnd.Position = [0.179464285714285,0.656547619047621,0.124107142857143,0.198809523809524];
    filename = [a_basePath filesep 'SynthImResp_allCells'];

else
    lgnd.Position = [0.20625,0.782142857142859,0.109821428571428,0.09047619047619];
    filename = [a_basePath filesep 'SynthImResp_highCorrCells'];

end
% print(f1, filename, '-dpng', '-r0')

%% normalize FRs and average - errorbar plot for all cells

cResp = [];

for idx = 1:length(aft_matches)
    
    cellId = cellfun(@(x) isequal(x, aft_matches(idx)), synthResponses(:, 3), 'UniformOutput', false);
    
    cellName = cell2mat(synthResponses(cell2mat(cellId), 3));
    
    cellResp = synthResponses{cell2mat(cellId), 1};
    cResp(:, idx) = cellResp((aft_idx(idx)-1)*num_steps+1:aft_idx(idx)*num_steps, 1);
    
    cResp(:, idx) = (cResp(:, idx) - mean(cResp(:, idx)))/std(cResp(:, idx));% - mean(cResp(:, idx)));
end

cResp = cResp';

err = std(cResp)./size(cResp, 1);

f2 = figure; 
hold on
title('Response to synthetic images along axis')
xlabel('Distance along preferred axis')
ylabel('Standardized firing rate')
e = errorbar(1:num_steps, mean(cResp, 1), err, 'k');
e.Marker = 'o';
e.MarkerFaceColor = 'k';
e.MarkerEdgeColor = 'k';
e.LineWidth = 2;
xtiklabs = {setdiff(-num_steps/2:1:num_steps/2, 0)};
xticks([1:20])
xticklabels(xtiklabs)
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

if ~exist('cells2use', 'var')
    filename = [a_basePath filesep 'SynthImResp_mean_allCells'];
else
    filename = [a_basePath filesep 'SynthImResp_mean_highCorrCells'];
end
print(f2, filename, '-dpng', '-r0')



%% error bar for single exemplar - using raster

num_trials = num_steps*4;


for idx = 1:length(aft_matches) % cell 1943
    
    cellId = cellfun(@(x) isequal(x, aft_matches(idx)), synthResponses(:, 3), 'UniformOutput', false);
    
    cellIndex = find(cell2mat(cellId) == 1);
    
    cellName = cell2mat(synthResponses(cell2mat(cellId), 3));
    
    
    windowBegin = synthResponses{cellIndex, 2} + 170;
    windowEnd = windowBegin+267;
    
    cellPsth = synthPsth{cellIndex, 1}(:, windowBegin:windowEnd);
    
    
    
    Resp = mean(cellPsth, 2)*1e3;
    cellResp = Resp((aft_idx(idx)-1)*num_trials+1:aft_idx(idx)*num_trials, 1);
    
    cellResp = reshape(cellResp, [4 20]);
    f3 = figure;
    hold on
    
    err = std(cellResp)./size(cellResp, 1);
    e = errorbar(1:num_steps, mean(cellResp, 1), err, 'LineWidth', 2);
    e.Color = cols(idx, :);
    e.Marker = 'o';
    e.MarkerFaceColor = cols(idx, :);
    e.MarkerEdgeColor = cols(idx, :);
    
    title({'Response to synthetic images', ['Axis corr ' num2str(axcorr(aft_idx(idx)))]})
    xtiklabs = {setdiff(-num_steps/2:1:num_steps/2, 0)};
    xticks([1:20])
    xticklabels(xtiklabs)
    lgnd = legend(legNames(aft_idx(idx)), 'FontWeight', 'bold');
    lgnd.Position = [0.20625,0.782142857142859,0.109821428571428,0.09047619047619];
    xlabel('Steps along preferred axis')
    ylabel('Spikes/s')
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
 
    if ~exist('cells2use', 'var')
        filename = [a_basePath filesep ['SynthImResp_' num2str(cellName)]];       
    else
        filename = [a_basePath filesep ['SynthImResp_exemplar_' num2str(cellName)]];        
    end
    print(f3, filename, '-dpng', '-r0')

    
end


close all

