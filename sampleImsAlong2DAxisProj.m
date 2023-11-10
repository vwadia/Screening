
% Script to automatically sample images along the 2D projection of preferred and orthogonal axes
% 
% for each cell
%     compute axis and orthogonal ax
%     get 2D scatter of projections
%     grab points along the beginning, middle, end of the xvalues
%     similar for y values
% vwadia Aug2022

%load stuff
setDiskPaths

taskPath = 'Object_Screening';
load([diskPath filesep taskPath filesep 'ITCells_500stim_Scrn_SigRamp']); 

options.screenType = 'Object';
stimDir = '500Stimuli';
layermat = 'fc6';
imageIDs = [1:500]';

load([diskPath filesep 'ObjectSpace' filesep stimDir filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu
params_full = scores{1}; 


[coeff, score, ~, ~, ~, mu] = pca(params_full);%, 'Centered', false);
params = score(:, 1:50);
    
options.screenType = 'Object';
options.pathStimuli = [diskPath filesep taskPath filesep '500Stimuli'];
options.ind_train = imageIDs;

pathOut = [diskPath filesep taskPath filesep 'SigRamp' filesep 'ImagesAlongAxes_2DProj'];
if ~exist(pathOut)
    mkdir(pathOut);
end
%%

for cellIndex = 1:length(strctCells)
    
    hfig = Utilities.ObjectSpace.sample_Ims_Along_STA(responses{cellIndex, 1}, params, options);
    sgtitle({['Cell number ' num2str(strctCells(cellIndex).Name)] 'Images along axes '});
    filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_' num2str(strctCells(cellIndex).Name) ];
    print(hfig, filename, '-dpng', '-r0')
    close all
    
end


%%

% 
% faceInds = 134:210;
% objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
% textInds = [264:282 290 400:408];
% vegInds = [211:235 357:399];
% animInds = [1:84 256:263];
% 
% 
% faceParams = params(faceInds, 1:2);
% textParams = params(textInds, 1:2);
% 
% 
% f = figure;
% hold on
% scatter(params(:, 1), params(:, 2), 50, [0 0 1], 'filled',  'HandleVisibility','off');
% scatter(faceParams(:, 1), faceParams(:, 2), 50, [1 0 0], 'filled');
% scatter(textParams(:, 1), textParams(:, 2), 50, [0 1 0], 'filled');
% xline(0);
% yline(0);
% % lg = legend({'Face Stimuli', 'Text Stimuli'});
% box off
% set(gca, 'Visible', 'off')
% print(f, 'PC12_FaceTextparamsSpread', '-dpng', '-r0')
