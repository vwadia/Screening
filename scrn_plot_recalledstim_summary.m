%% Plotting the encoding spikerate of recalled stimuli
%% raster + bargraph side by side

% Take total raster and order
% find stim Raster
% Plot it in one subplot
% Make a barplot for mean spikes/s in the appropriate window on the other
% create a figure that represents the order the stimuli appear on the selected cells STA
% Needs strctCells, & screeningData

%% 
color = [0 0.2 0];
MarkerSize = 4;
spikeCount = [];

% collect images
imDir = dir(fullfile(pathStimuli));
imDir = imDir(~ismember({imDir.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));


% options.recalled_stim = [12 19 25 123 270 487]; % P76 Recall 1
options.recalled_stim = [54 129 130 186 270 449]; % P76 Recall 2
for cellIndex = l(strctCells)
     if ~isempty(screeningData.responses{cellIndex, 2}) && isfield(options, 'recalled_stim')
         for rs = 1:length(options.recalled_stim)
             stimRaster = screeningData.psth{cellIndex, 1}(find(screeningData.sortedOrder == options.recalled_stim(rs)), :);
             offsetSC = screeningData.responses{cellIndex, 2}-screeningData.timelimits(1)*1e3;
             
             if size(stimRaster, 2) >= offsetSC
                 spikeCount(rs, 1) = mean(mean(stimRaster(:, offsetSC:offsetSC+ceil(stimDur))))*1e3; % mulitply by 1e3 to make it spikes/s                                  
             else
                 spikeCount(rs, 1) = mean(mean(stimRaster(:, offsetSC:offsetSC+ceil(stimDur))))*1e3; % mulitply by 1e3 to make it spikes/s                 
             end
         end 
         
         
         for rs = 1:length(options.recalled_stim)
             % compute stimraster
             stimRaster = screeningData.psth{cellIndex, 1}(find(screeningData.sortedOrder == options.recalled_stim(rs)), :);
             
             % relevant variables
             imPath = [imDir(options.recalled_stim(rs)).folder filesep imDir(options.recalled_stim(rs)).name];
             pathOut = [basePath filesep 'rasters' filesep 'EncodingSpikeCount'];
             if ~exist(pathOut)
                 mkdir(pathOut);
             end
             filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber)...
                 '_' num2str(strctCells(cellIndex).Name) '_Stim_' num2str(rs)];
             plotOptions.timelimits = screeningData.timelimits;
             plotOptions.globalyl = max(spikeCount);
             
             % produce figure 
             handlesToFig = Utilities.Plotting.PlotRasterAndBar_SingleStim(stimRaster, spikeCount(rs, 1), imPath, plotOptions);
             sgtitle({[strctCells(cellIndex).brainArea ' ' num2str(strctCells(cellIndex).Name)],...
                 ['Encoding Spike Count Stim ' num2str(rs)]})
             
             % save
             print(handlesToFig, filename, '-dpng', '-r0');
             close all
         end
     end
end

