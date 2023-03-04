setDiskPaths

% load([diskPath filesep 'Object_Screening' filesep 'AllRespITCells_Morn&AftSessions_withPDist_Scrn_500Stim.mat'])
load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat'])

%% orange plot all cells - categories arranged as in legend of rasters
imageIDs = [1:500]';

faceInds = 134:210;
objInds = [85:133 236:255 283:356 409:500];
textInds = [264:282 400:408];
vegInds = [211:235 357:399];
animInds = [1:84 256:263];

allInds = zeros(length(imageIDs), 1);
allInds(ismember(imageIDs, faceInds)) = 1;
allInds(ismember(imageIDs, textInds)) = 2;
allInds(ismember(imageIDs, vegInds)) = 3;
allInds(ismember(imageIDs, animInds)) = 4;
allInds(ismember(imageIDs, objInds)) = 5;

[s, c] = sortrows(allInds);

blineStim = NaN(length(strctCells), 1);
respMat = [];
for cellIndex = 1:length(strctCells)
    if ~isempty(responses{cellIndex, 2})
        
        % calculate baseline per cell
        respLat = responses{cellIndex, 2};
        if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
            sortedOrder = repelem(imageIDs, 6);
            stimDur = 166.6250;
            stimOffDur = 166.6250;
            timelimits = [-0.17, 0.33];
        else
            sortedOrder = repelem(imageIDs, 4);
            stimDur = 266.6250;
            stimOffDur = 133.500; 
            timelimits = [-0.13, 0.53];
        end
        windowBegin = (-timelimits(1)*1e3)+responses{cellIndex, 2};    
        if stimOffDur > windowBegin
            blineStim(cellIndex, 1) = mean(mean(psths{cellIndex, 1}(:, 1:windowBegin)));
        else
            blineStim(cellIndex, 1) = mean(mean(psths{cellIndex, 1}(:, windowBegin-ceil(stimOffDur):windowBegin)))*1e3;
        end
        respMat = [respMat responses{cellIndex, 1}(c, 1)];
        
    end
end
% blineStim = blineStim(~isnan(blineStim));
respMat = respMat';

% now use calculated values to Normalize
respMatNorm = respMat - repmat(blineStim, 1, size(respMat, 2));

for cellIndex = 1:size(respMat, 1)
    mx = max(abs(respMatNorm(cellIndex, :)));
    respMatNorm(cellIndex, :) = respMatNorm(cellIndex, :)./mx;
end

% respMatNorm = respMatNorm(:, c);
%% sorting cells by max category

for cellIndex = 1:length(strctCells)
    
    % sum the max category response
    for cat = 1:5
        catr(cellIndex, cat) = sum(respMatNorm(cellIndex, s == cat));
    end
    [~, maxCat(cellIndex, 1)] = max(catr(cellIndex, :));
    
end
 
[a, b] = sort(maxCat);
respMatNorm = respMatNorm(b, :);


%% Plot orangemap
cols = Utilities.distinguishable_colors(5);
cols = cols([2 3 4 5 1], :)
f = figure; clf;
set(gcf,'Position',get(0,'Screensize'))
% hold on
I = imagesc(respMatNorm, [-max(max(respMatNorm)), max(max(respMatNorm))]);

ylabel('Cell Number', 'FontSize', 16, 'FontWeight', 'bold');
ytik = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', ytik, 'fontsize', 16, 'FontWeight', 'bold');

tit = title([{'Population summary of all responsive units - IT cortex'}]);
tit.FontSize = 20; 

% setting up my colormap
orangemap = esa(300);
[WhiteColor, Whitepos] = max(sum(orangemap,2));
orangemap = orangemap([1:Whitepos round(linspace(Whitepos+1,size(orangemap,1)-2,Whitepos-1))],:);
colormap(orangemap)

colbar = colorbar;
colbar.AxisLocation = 'out';
colbar.Position = [0.935119047619047,0.702380952380952,0.020238095238096,0.21666666666667];
colbar.Ticks = [-max(max(respMatNorm)), max(max(respMatNorm))];
colbar.TickLabels = {'-1', '1'};
colbar.Label.String = 'Norm firing rate';

% setting xticklabels correctly
xtikpos = zeros(1, length(unique(allInds)));
ctr = 1;
for i = unique(allInds)'
    if i == 1
        xtikpos(ctr) = length(find(allInds == i));
    else
        xtikpos(ctr) = xtikpos(ctr-1) + length(find(allInds == i));
    end
    ctr = ctr + 1;
end
xticks(xtikpos);
set(gca, 'xticklabel', {[]}); % clear the old labels
xtik = get(gca,'xtick'); 
xtiklabs = {'Faces', 'Text', 'Plants/fruits', 'Animals', 'Objects'};

for ticknum = 1:length(xtik)
   if ticknum == 1
       xtikpos = xtik(ticknum)/3;
   else
       xtikpos=xtik(ticknum-1)+(xtik(ticknum)-xtik(ticknum-1))/4;
   end
   tx = text(xtikpos, size(respMatNorm, 1)+5, xtiklabs{ticknum});
   tx.FontWeight = 'bold';
   tx.FontSize = 16;
%    set(tx, 'Rotation', 45);
end
% marking categories
for i = 1:5
    rectangle('Position',[length(find(allInds == 0:i-1))+0.5 length(find(a == 0:i-1))+0.5 length(find(allInds == i))  length(find(a == i))], 'EdgeColor',cols(i, :), 'LineWidth', 3);
end
filename = [diskPath filesep 'Object_Screening' filesep 'OrangePlot_AllRespUnits_IT_markedCats'];
print(f, filename, '-dpng', '-r0')

%% RESPONSE LATENCY
%% histogram of response latency

% load([diskPath filesep 'Object_Screening' filesep 'ITRespCells_500stim_Scrn.mat']); 
respLat = cell2mat(responses(:, 2));
f = figure;
hold on
histogram(respLat, 50:25:375, 'FaceColor', [0.6350 0.0780 0.1840], 'EdgeColor', [0 0 0]); 
filename = [diskPath filesep 'Object_Screening' filesep 'Hist_RespLat_IT'];
title('Response latency of IT neurons', 'FontSize', 16, 'FontWeight', 'bold')


% histogram(respAmg, 25:50:650, 'FaceColor', [0.4940 0.1840 0.5560], 'EdgeColor', [0 0 0]); 
% filename = ['G:\SUAnalysis\Object_Screening\Hist_RespLat_ITvsAmyg'];
% title('Response latency of IT vs Amyg neurons', 'FontSize', 16, 'FontWeight', 'bold')



xticks([50:50:600]);
xlabel('Response Latency (ms)')
ylabel('Number of cells')
set(gca, 'FontWeight', 'bold');
% ylim([0 40])
% lgnd = legend({'IT', 'Amyg'}, 'FontWeight', 'bold');
% legend
print(f, filename, '-dpng', '-r0')

%% response latency comparison

% responses = screeningData
for i = 1:length(responses)
    
    if isequal(strctCells(i).Name, cell2mat(responses(i, 3)))
        responses{i, 4} = strctCells(i).brainArea;
    end
    
end

respLats = {};

ars = unique(responses(:, 4));
for ar = 1:length(ars)
    derp = cellfun(@(x) strcmp(ars{ar}, x), responses(:, 4), 'UniformOutput', false);
    
    respLats{ar, 1} = responses(cell2mat(derp), 2);
    
end
respLats(:, 2) = ars;
respAmg = cell2mat([respLats{1, 1}; respLats{5, 1}]);
% for i = 1:length(responses)
%     if strcmp(responses{i, 4}, 'LSMA') || strcmp(responses{i, 4}, 'RSMA')
%         respLat{i, 1} = responses{i, 2};
%     elseif strcmp(responses{i, 4}, 'LA') || strcmp(responses{i, 4}, 'RA')
%     elseif strcmp(responses{i, 4}, 'LH') || strcmp(responses{i, 4}, 'RH')
%     elseif strcmp(responses{i, 4}, 'LAC') || strcmp(responses{i, 4}, 'RAC')
%     elseif strcmp(responses{i, 4}, 'LOF') || strcmp(responses{i, 4}, 'ROF')
%     elseif strcmp(responses{i, 4}, 'RFFA') || strcmp(responses{i, 4}, 'LFFA')
%     end
% end


