% Viewpoint plot per identity

respMatNormFV = respMatNorm(:, 1:200);
%% Per identity
pathOut = [basePath filesep 'rasters' filesep 'ViewPoint'];
if ~exist(pathOut)
    mkdir([pathOut]);
end
cols = distinguishable_colors(8);
for cellIndex = 1:size(respMatNormFV, 1) % cycle through cells
    vec = respMatNormFV(cellIndex, :);
    f = figure; clf;
    set(gcf,'Position',get(0,'Screensize'))
    hold on
    
    for idx = 1:8
        plot(vec(find(s == idx)), '-o', 'Color', cols(idx, :));
    end
    
    lgnd = legend('front', 'halfLeft', 'fullLeft', 'halfRight', 'fullRight', 'lookUp', 'lookDown', 'backOfHead');
    title([num2str(strctCells(cellIndex).ChannelNumber) '\_' num2str(strctCells(cellIndex).Name) '\_' char(strctCells(cellIndex).brainArea) '\_ViewPoint']);
    ylabel('Norm Firing rate', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Identity #', 'FontSize', 16, 'FontWeight', 'bold');
    xticks([1:25]);
    set(get(gca, 'XAxis'), 'FontWeight', 'bold');
    set(get(gca, 'YAxis'), 'FontWeight', 'bold');
    xlim([1 25]);
    ylim([-1 1]);
    set(lgnd, 'Position', [0.919,0.770,0.061,0.145]);

    filename = [pathOut filesep strctCells(cellIndex).brainArea '_' num2str(strctCells(cellIndex).ChannelNumber) '_'...
        num2str(strctCells(cellIndex).Name) '_PerIdentity'];
    print(f, filename, '-dpng', '-r0')
    close all
end

%% make RSA matrices for viewpoints - need to collect viewpoints properly

viewMat = [];

for cellIndex = l(strctCells)
    if ~isempty(screeningData.responses{cellIndex, 2})
        viewMat = [viewMat screeningData.responses{cellIndex, 1}]; 
    end
end

viewMat = viewMat(1:200, :); % get rid of non-face objects

% collect all 25 of each view together
sm_order = [];
new_view = [];
for i = 1:8
    new_view  = [i:8:200];
    if i == 1
        sm_order = new_view;
    else
        sm_order = [sm_order new_view];
    end
end

viewMat = viewMat(sm_order, :); % arrange matrix accordingly
viewMat = viewMat';
corrMat = corr(viewMat);
%%
f = figure; clf;
% set(gcf,'Position',get(0,'Screensize'))
hold on
colormap(flipud(gray));
imagesc(corrMat);
colbar = colorbar;
colbar.Label.String = 'Correlation';
colbar.FontWeight = 'Bold';
set(get(gca, 'XAxis'), 'Fontweight', 'Bold');
set(get(gca, 'YAxis'), 'Fontweight', 'Bold');

xlim([0 200])
ylim([0 200])
xticks([0 25 50 75 100 125 150 175 200])
yticks([0 25 50 75 100 125 150 175 200])

filename = [basePath filesep 'SimilarityMatrixViewpoints'];
print(f, filename, '-dpng', '-r0');
close all


