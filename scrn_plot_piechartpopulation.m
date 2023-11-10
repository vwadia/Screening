% Population plotting for FFA screens
% This script needs
%     screeningData struct with (raster, psth, times)
%     strctCells
%     Sort the images by category 
%     Sort the responses by Face Selectivity index (FSI)
%     make orange plot

% Note: when subtracting baseline and dividing by max one can do this on a
% per stim bases. I did not do that here, kept 1 baseline value for the
% cell across trials (See face procesing notes on iPad)

% vwadia April2021

setDiskPaths

load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim'])


%%
% get rid of P71's cells 
% only for extraction of face cells
index = zeros(length(strctCells), 1);

for id = 1:length(strctCells)
    
    if strcmp(strctCells(id).PatientID, 'P71CS')
        index(id) = 1;
    end
    
end

index = logical(index);
responses(index, :) = [];
psths(index, :) = [];
strctCells(index) = [];


screeningData = struct;
screeningData.imageIDs = [1:500]';
screeningData.psth = psths;
screeningData.responses = responses;
screeningData.sortedOrder = repelem([1:500], 4)';

%% 500 object screen

faceInds = 134:210;
objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
textInds = [264:282 290 400:408];
vegInds = [211:235 357:399];
animInds = [1:84 256:263];


allInds = zeros(length(screeningData.imageIDs), 1);
allInds(ismember(screeningData.imageIDs, faceInds)) = 1;
allInds(ismember(screeningData.imageIDs, textInds)) = 2;
allInds(ismember(screeningData.imageIDs, vegInds)) = 3;
allInds(ismember(screeningData.imageIDs, animInds)) = 4;
allInds(ismember(screeningData.imageIDs, objInds)) = 5;

%% 1593 object screen

% faceInds = 412:612;
% objInds = [217:411 694:768 874:1122 1246:1593];
% textInds = [796:873 1220:1245];
% vegInds = [613:693 1123:1219];
% animInds = [1:216 769:795];
% 
% allInds = zeros(length(screeningData.imageIDs), 1);
% allInds(ismember(screeningData.imageIDs, faceInds)) = 1;
% allInds(ismember(screeningData.imageIDs, textInds)) = 2;
% allInds(ismember(screeningData.imageIDs, vegInds)) = 3;
% allInds(ismember(screeningData.imageIDs, animInds)) = 4;
% allInds(ismember(screeningData.imageIDs, objInds)) = 5;

%% Category screen - indexing to sort images by category

% allInds = zeros(length(screeningData.imageIDs), 1);
% allInds(ismember(screeningData.imageIDs, faceInds)) = 1;
% allInds(ismember(screeningData.imageIDs, textInds)) = 2;
% allInds(ismember(screeningData.imageIDs, placeInds)) = 3;
% allInds(ismember(screeningData.imageIDs, scrambledInds)) = 4;
% allInds(ismember(screeningData.imageIDs, objInds)) = 5;

%% Face view screen

% Separating faces and objects (combining views)
% allInds = zeros(length(screeningData.imageIDs), 1);
% allInds(ismember(screeningData.imageIDs, front)) = 1;
% allInds(ismember(screeningData.imageIDs, halfLeft)) = 1;
% allInds(ismember(screeningData.imageIDs, fullLeft)) = 1;
% allInds(ismember(screeningData.imageIDs, halfRight)) = 1;
% allInds(ismember(screeningData.imageIDs, fullRight)) = 1;
% allInds(ismember(screeningData.imageIDs, lookUp)) = 1;
% allInds(ismember(screeningData.imageIDs, lookDown)) = 1;
% allInds(ismember(screeningData.imageIDs, backOfHead)) = 9;
% allInds(ismember(screeningData.imageIDs, objects)) = 9;

% separating views
% allInds = zeros(length(screeningData.imageIDs), 1);
% allInds(ismember(screeningData.imageIDs, front)) = 1;
% allInds(ismember(screeningData.imageIDs, halfLeft)) = 2;
% allInds(ismember(screeningData.imageIDs, fullLeft)) = 3;
% allInds(ismember(screeningData.imageIDs, halfRight)) = 4;
% allInds(ismember(screeningData.imageIDs, fullRight)) = 5;
% allInds(ismember(screeningData.imageIDs, lookUp)) = 6;
% allInds(ismember(screeningData.imageIDs, lookDown)) = 7;
% allInds(ismember(screeningData.imageIDs, backOfHead)) = 8;
% allInds(ismember(screeningData.imageIDs, objects)) = 9;

% separating identities
% allInds = zeros(length(screeningData.imageIDs), 1);
% allInds = [1:25];
% allInds = repelem(allInds, 8)';
% allInds = [allInds; repelem(99, 32)'];

[s, c] = sortrows(allInds);
%% calculate responses to each image - not that these are arranged so that faces (category 1) are first - hence different from screeningData.responses 

% use individual image information
orderToUse = screeningData.sortedOrder;
totalTrials = length(screeningData.imageIDs);
sortedLabels = c; % screeningData.imageIDs;

blineStim = zeros(length(strctCells), 1);
respMat = NaN(length(strctCells), length(sortedLabels));

offset = 0;

for cellIndex = 1:length(strctCells)
    
if strcmp(strctCells(cellIndex).SessionID, 'P79CS_2')...
            || strcmp(strctCells(cellIndex).SessionID, 'P79CS_3')...
            || strcmp(strctCells(cellIndex).SessionID, 'P79CS_4')
%         orderToUse = repelem([1:500], 4)';
        stimDur = 266.6250;
        windowLength = ceil(stimDur);
        screeningData.timelimits = [-0.17 0.53];
        stimOffDur = 166.6250;
    else
%         orderToUse = repelem([1:500], 4)';
        stimDur = 266.6250;
        windowLength = ceil(stimDur);
        screeningData.timelimits = [-0.13 0.53];
        stimOffDur = 166.6250;
    end
    
    if ~isempty(screeningData.responses{cellIndex, 2})
        windowBegin = (-screeningData.timelimits(1)*1e3)+screeningData.responses{cellIndex, 2};     
        
        % single number for baseline per cell
        % WTF IS THIS EVEN DOING???
        if stimOffDur > windowBegin
            blineStim(cellIndex, 1) = mean(mean(screeningData.psth{cellIndex, 1}(:, 1:windowBegin)));
        else
            blineStim(cellIndex, 1) = mean(mean(screeningData.psth{cellIndex, 1}(:, windowBegin-ceil(stimOffDur):windowBegin)));
        end
        
        
        for img = 1:length(sortedLabels)
            stimMat = screeningData.psth{cellIndex, 1}(find(orderToUse == sortedLabels(img)), :);
            % response of a cell to a given stimulus
            if windowBegin+windowLength > size(stimMat, 2)
                respStim = mean(mean(stimMat(:, windowBegin:end)));
            else
                respStim = mean(mean(stimMat(:, windowBegin:windowBegin+windowLength)));
            end
            respMat(cellIndex, img) = respStim;
        end
    end
end



% now use calculated values to Normalize
respMatNorm = respMat - repmat(blineStim, 1, size(respMat, 2));

for cellIndex = l(strctCells)
    mx = max(abs(respMatNorm(cellIndex, :)));
    respMatNorm(cellIndex, :) = respMatNorm(cellIndex, :)./mx;
end

%% Calculate FSI

% Face Selectivity Index - per cell
FSI = NaN(length(strctCells), 1);
faceResp = 0;
nonFaceResp = 0;
faceIdx = 1; % set in first cell

[screeningData.catOrder, ~] = Utilities.makeObjCatLabsScreening('P76CSFast', orderToUse);      



for cellIndex = l(strctCells)
%     if cellIndex == 22
%         keyboard
%     end
%     if strcmp(strctCells(cellIndex).PatientID, 'P71CS')
%         orderToUse = repelem([1:500], 6)';
%         screeningData.catOrder = Utilities.makeObjCatLabsScreening(strctCells(cellIndex).SessionID, orderToUse);
%         stimDur = 166.6250;
%         windowLength = ceil(stimDur);
%         screeningData.timelimits = [-0.17 0.33];
%         stimOffDur = 166.6250;
%     if strcmp(strctCells(cellIndex).SessionID, 'P79CS_2')...
%             || strcmp(strctCells(cellIndex).SessionID, 'P79CS_3')...
%             || strcmp(strctCells(cellIndex).SessionID, 'P79CS_4')
%         orderToUse = repelem([1:500], 4)';
%         screeningData.catOrder = Utilities.makeObjCatLabsScreening(strctCells(cellIndex).SessionID, orderToUse);       
%         screeningData.timelimits = [-0.17 0.53];
%     else
%         orderToUse = repelem([1:500], 4)';
%         screeningData.catOrder = Utilities.makeObjCatLabsScreening(strctCells(cellIndex).SessionID, orderToUse);      
%         screeningData.timelimits = [-0.13 0.53];
%     end
    
    
    if ~isempty(screeningData.responses{cellIndex, 2})
        windowBegin = (-screeningData.timelimits(1)*1e3)+screeningData.responses{cellIndex, 2};
        
        faceMat = screeningData.psth{cellIndex, 1}(find(screeningData.catOrder == faceIdx), :);
        nonFaceMat = screeningData.psth{cellIndex, 1}(find(screeningData.catOrder ~= faceIdx), :);
        
        if windowBegin+windowLength > size(faceMat, 2)
            faceResp = mean(mean(faceMat(:, windowBegin:end), 1));
            nonFaceResp = mean(mean(nonFaceMat(:, windowBegin:end), 1));
        else
            faceResp = mean(mean(faceMat(:, windowBegin:windowBegin+windowLength), 1));
            nonFaceResp = mean(mean(nonFaceMat(:, windowBegin:windowBegin+windowLength), 1));
        end
        FSI(cellIndex) = (faceResp - nonFaceResp)/(faceResp + nonFaceResp);
    end
end

faceCellInds = find(FSI > 0.33);

%% save face cells - for ev calculation

responses = screeningData.responses;
psths = screeningData.psth;
sortedOrder = screeningData.sortedOrder;
catOrder = screeningData.catOrder;

responses = responses(faceCellInds, :);
psths = psths(faceCellInds, :);
faceCells = strctCells(faceCellInds);

% save([diskPath filesep 'Object_Screening' filesep 'AllFaceCells_1593_Objects'], 'faceCells', 'responses', 'psths', 'catOrder', 'sortedOrder', '-v7.3');
save([diskPath filesep 'Object_Screening' filesep 'FaceCells_500Stim_Scrn'], 'faceCells', 'responses', 'psths', 'catOrder', 'sortedOrder', '-v7.3');

%%
[s_FSI, c_FSI] = sortrows(FSI, 'descend');

screeningData.strctCells_sorted_FSI = strctCells(c_FSI); 
screeningData.strctCells_sorted_FSI = screeningData.strctCells_sorted_FSI(~isnan(s_FSI));
%% distribution of FSI

f = figure; clf;
hold on
hist(FSI);
xline(1/3, '--', 'LineWidth', 2);
xline(-1/3, '--', 'LineWidth', 2);
xlim([-1 1]);
ylim([0 10]);
box on
xlabel('Face Selectivity Index (FSI)')
ylabel('Number of cells')
title('Distribution of FSI Face View Screen March28');
filename = [basePath filesep 'FSIValueDist_2'];
print(f, filename, '-dpng', '-r300')
close all
%%

f = figure; clf;
set(gcf,'Position',get(0,'Screensize'))
hold on

% flipud because I want the cells to be plotted from the topdown (imageSC
% plots from the bottom up)
I = imagesc(respMatNorm(flipud(c_FSI(~isnan(s_FSI))), :), [-max(max(respMatNorm)), max(max(respMatNorm))]);

tit = title([{'Category screening'}]);
tit.FontSize = 20; 

ylabel('Cell Number', 'FontSize', 16, 'FontWeight', 'bold');
ytik = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', ytik, 'fontsize', 16, 'FontWeight', 'bold');
% xlabel('Images');
ylim([0 size(c_FSI(~isnan(s_FSI)), 1)+1]);
yl = ylim;
xlim([0 length(sortedLabels)+1]);

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
xtiklabs = screeningData.lgnd;

for ticknum = 1:length(xtik)
   if ticknum == 1
       xtikpos = xtik(ticknum)/3;
   else
       xtikpos=xtik(ticknum-1)+(xtik(ticknum)-xtik(ticknum-1))/4;
   end
   tx = text(xtikpos, -1, xtiklabs{ticknum});
   tx.FontWeight = 'bold';
   tx.FontSize = 16;
%    set(tx, 'Rotation', 45);
end


destPath = basePath;
filename = [destPath filesep [taskStruct.filePrefix '_sortedFSI']];

print(f,filename ,'-dpng','-r300')
close all
%% pie chart - inputted manually

% Oct 2022 
% total cell count = 334 (284 vis resp/410 total - after merge have 235 resp cells - so scaled up the number to keep percentage intact as can't merge the non-resp cells 
% resp cells - 235
% non_resp cells - 99
% ramp cells - 165
% non ramp cells - 70
setDiskPaths

Total_merged_cells = 440; % note this is a guestimate
% load all resp neurons
load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim'])
respCells = length(strctCells);

load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500Stim_Scrn_SigRamp'])
sigRampCells = length(strctCells);


% create vectors of categorical labels
pieData = [sigRampCells respCells-sigRampCells Total_merged_cells-respCells];
% cats = {};
% for i = 1:length(pieData)
%     cats = horzcat(cats, repmat(labels(i), [1 pieData(i)]));
% end
%%
setDiskPaths
fig = figure; 
% set(gcf,'Position',get(0,'Screensize'))

% colormap bone
explode = [1 1];
% colormap turbo
% explode = [1 1 1 1 1];

% labels = {'Faces (26.04%)','Text (5.92%)', 'Plants/fruits (5.33%)' ,'Animals (23.67%)', 'Objects (39.05%)'};
% txt1 = ['Sig Ramp Tuned'  '(' num2str((sigRampCells/Total_merged_cells)*1e2,' %.2f') '%)'];
% txt2 = ['Responsive only' '(' num2str(((respCells-sigRampCells)/Total_merged_cells)*1e2,' %.2f') '%)'];
% txt3 = ['Neither' '(' num2str(((Total_merged_cells-respCells)/Total_merged_cells)*1e2, ' %.2f') '%)'];
% labels = {txt1, txt2, txt3};

% for imagination hand inputted
% colormap turbo
% explode = [1 1];
% pieData = [1019-366 366];
% txt1 =['Non-active during imagery', ' (' num2str(((1019-366)/1019)*1e2,' %.2f') '%)'];
% txt2 = ['Active during imagery', ' (' num2str((366/1019)*1e2,' %.2f') '%)'];
% labels = {txt1, txt2};

pos(1, :) = [0.15864406779661,1.114576271186441,0];
pos(2, :) = [0.211229540688194,-1.167055515514338,0];
pos(3, :) = [0.15864406779661,1.114576271186441,0];
% pos(4, :) = [0.211229540688194,-1.167055515514338-1,0];

% pieData = [173-75 44 31];
% pieData = [451 54];
pieData = [46 127];

% 451 - resp neurons merged, 534 - resp neurons unmerged  
% 597 - total neurons unmerged, estimate 505 - total neurons merged  

% apple = pie(pieData,explode, labels) 
apple = pie(pieData,explode); 
% for i = 2:2:length(apple)
for i = 2:2:length(apple)
    txt = apple(i);
%     txt.Position = pos(i/2);
    txt.FontSize = 14;
    txt.FontWeight = 'Bold';
end
% apple(8).Position = [0.55 1.1 0]; % manual lebal adjustment

% print(fig, [diskPath filesep 'Object_Screening' filesep 'Screening_Resp_Cells_pie_chart'], '-dpng', '-r300')
print(fig, [diskPath filesep 'Recall_Task' filesep 'Im_React_Cells_pie_chart'], '-dpng', '-r300')
% print(fig, [diskPath filesep 'Object_Screening' filesep 'Ramp_cells_pie_chart'], '-dpng', '-r300')
% close all


% print(fig, [diskPath filesep 'Recall_Task' filesep 'forPaper' filesep 'pieChart_IT_reac'], '-dpng', '-r300')
% close all











