
% This script loads in all IT cells ever recorded
% loops through them, sets paths and grabs the correct order and events for the session
% makes the psths and then compute response latency and visual responsivity

%% TAKEAWAY
% There are many more sigRamp cells if one uses sel latency instead of resp latency
% ~121 vs ~105

% if you use fc7 post relu you get the most consistent answers
% 120 sig ramp (sel resp) vs 103-105 (vis resp)



%% setpaths and load in IT cells

setDiskPaths
taskPath = 'Object_Screening';

% need to do everything from the server - this address is the same for
% connection server from lab desktop and computer at Cedars
addpath(genpath(['Code' filesep 'osortTextUI']));
addpath(genpath('ObjectSpace'));
% addpath(genpath('synthetic_face_generator'));



%%
load([diskPath filesep 'Object_Screening' filesep 'AllITCells_500stim_Scrn.mat'])

m_psths = psths;
m_responses = responses;
m_strctCells = strctCells;

load([diskPath filesep 'Object_Screening' filesep 'AllITCells_500stim_ReScreen.mat'])

a_psths = psths;
a_responses = responses;
a_strctCells = strctCells;

strctCells = [m_strctCells a_strctCells];
psths = [m_psths; a_psths];
responses = [m_responses(:, 1:3); a_responses(:, 1:3)];

ITCellsOnly = 1;

strctCELL = struct2cell(strctCells');
strctCELL = strctCELL';


%% set up orders

stimDir = '500Stimuli';
imageIDs = [1:500]';
order1 = repelem(imageIDs, 6);
order2 = repelem(imageIDs, 4);
catOrd = {};

% afternoon session that got fucked
a_basePath = [diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917'];
a_strct = load([a_basePath filesep 'PsthandResponses']);
order3 = a_strct.order;
clearvars a_strct


% make cat labels
anovaType = 'CategoryObject';
faceInds = 134:210;
objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
textInds = [264:282 290 400:408];
vegInds = [211:235 357:399];
animInds = [1:84 256:263];

for i = 1:3
    
    if i == 1
        order = order1; % 3000 elements
    elseif i == 2
        order = order2; % 2000 elements
    elseif i == 3
        order = order3;
    end
    
    catOrder = zeros(length(order), 1);
    catOrder(ismember(order, faceInds)) = 1;
    catOrder(ismember(order, textInds)) = 2;
    catOrder(ismember(order, vegInds)) = 3;
    catOrder(ismember(order, animInds)) = 4;
    catOrder(ismember(order, objInds)) = 5;
    
    catOrd{i} = catOrder;
end
%% P73CS FullParam_Obj Screen
%  load([diskPath filesep taskPath filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328'...
%         filesep 'ITCells_1593Stim_Scrn'])
%
%  load([diskPath filesep taskPath filesep 'P73CS' filesep 'FullParamObjScreening_Session_1_20210328'...
%     filesep 'realOrder_P73CS_ParamObj']); % produces imageTexture
%
% stimDir = '1593Stimuli';
% order = imageTexture - 10; order = order';
%
% faceInds = 412:612;
% objInds = [217:411 694:768 874:1122 1246:1593];
% textInds = [796:873 1220:1245];
% vegInds = [613:693 1123:1219];
% animInds = [1:216 769:795];
%
% catOrder = zeros(length(order), 1);
% catOrder(ismember(order, faceInds)) = 1;
% catOrder(ismember(order, textInds)) = 2;
% catOrder(ismember(order, vegInds)) = 3;
% catOrder(ismember(order, animInds)) = 4;
% catOrder(ismember(order, objInds)) = 5;
% catOrd = catOrder;
%
% imageIDs = [1:1593]';
%% compute new responses with different ways and compare

% which to use
% 0 - basic,
% 1 - sliding window anova,
% 2 - peak of omega squared (also sliding window anova)
% 3 - Poisson single trial latency (added March 2023 vwadia)

% March2023 -
% 0 - 365/226 responsive/ramp tuned units respectively
% 1 - 237/207 responsive/ramp tuned units respectively
% 2 - 329/201 responsive/ramp tuned units respectively
% 3 - 397/240 responsive/ramp tuned units respectively

c_resp = {};

comparison(:, 1) = strctCELL(:, 1);

for method = 3%[0 1 2 3]
    % method = 0;
    compResponses = cell(length(strctCells), 3);
    tic
    for cellIndex = l(strctCells)
        
        if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
            labels = catOrd{1};
            timelimits = [-0.17, 0.33];
            stimOffDur = 166.6250;
            stimDur = 166.6250;
            sortedOrder = order1;
        elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
            labels = catOrd{3};
            timelimits = [-0.17, 0.33];
            stimOffDur = 133.4680;
            stimDur = 266.6250;
            sortedOrder = order3;
        else
            labels = catOrd{2};
            timelimits = [-0.17, 0.53]; %changed for all cells Oct2022 vwadia
            stimOffDur = 133.4680;
            stimDur = 266.6250;
            sortedOrder = order2;
        end
        
        %      labels = catOrd;
        %      timelimits = [-0.13, 0.53];
        %      stimOffDur = 133.4680;
        %      stimDur = 266.6250;
        %      [sortedOrder, correctOrder] = sortrows(order);
        if method ~= 3
            n_stdDevs = 2.5;
            [respLat, max_group] = Utilities.computeResponseLatency(psths(cellIndex, :), labels, timelimits,...
                stimOffDur, stimDur, method, n_stdDevs);
            if isempty(respLat)
                respLat = 0;
            end
            endRas = size(psths{cellIndex, 1}, 2);
            if respLat ~= 0
                
                windowLength = floor(stimDur);
                windowBegin = respLat;
                windowEnd = windowBegin+windowLength;
                if windowEnd > endRas
                    windowEnd = endRas;
                end
                
                
                for i = l(imageIDs)
                    stimRaster = psths{cellIndex, 1}(find(sortedOrder == imageIDs(i)), windowBegin:windowEnd);
                    compResponses{cellIndex, 1}(i, 1) = mean(mean(stimRaster))*1e3; % note the 1000x multiplcation
                end
                
                compResponses{cellIndex, 2} = respLat - (-timelimits(1)*1e3);
                compResponses{cellIndex, 3} = strctCells(cellIndex).Name;
                
                %         compResponses{cellIndex, 4} = strctCells(cellIndex).ChannelNumber;
                %         compResponses{cellIndex, 5} = strctCells(cellIndex).SessionID;
                
                if exist('max_group')
                    compResponses{cellIndex, 4} = max_group;
                end
                
                if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
                    if compResponses{cellIndex, 2} > 250
                        compResponses(cellIndex, :) = [];
                    end
                elseif ITCellsOnly && compResponses{cellIndex, 2} > 350
                    compResponses(cellIndex, :) = [];
                end
                
                respLat = 0; % empty it so it gets assigned again if present
                
            end
        else
            [respLat, ~] = Utilities.computeRespLatPoisson(psths(cellIndex, :), labels, sortedOrder, timelimits,...
                stimDur, true); % 3.5 std devs and trial baseline
            endRas = size(psths{cellIndex, 1}, 2);
            
            if ~isnan(respLat)
                
                windowLength = floor(stimDur);
                windowBegin =  (-timelimits(1)*1e3) + floor(respLat);
                windowEnd = windowBegin+windowLength;
                if windowEnd > endRas
                    windowEnd = endRas;
                end
                
                
                for i = l(imageIDs)
                    stimRaster = psths{cellIndex, 1}(find(sortedOrder == imageIDs(i)), windowBegin:windowEnd);
                    compResponses{cellIndex, 1}(i, 1) = mean(mean(stimRaster))*1e3; % note the 1000x multiplcation
                end
                
                compResponses{cellIndex, 2} = respLat;
                compResponses{cellIndex, 3} = strctCells(cellIndex).Name;
                
                %         compResponses{cellIndex, 4} = strctCells(cellIndex).ChannelNumber;
                %         compResponses{cellIndex, 5} = strctCells(cellIndex).SessionID;
                
                if exist('max_group')
                    compResponses{cellIndex, 4} = max_group;
                end
                
            end
            
            
        end
        %     if exist('max_group')
        %         compResponses{cellIndex, 4} = max_group;
        %     end
    end
    toc
    
    rL{method+1, 1} = cell2mat(compResponses(:, 2));
    disp(length(rL{method+1, 1}))
    
    % checking if pruned cells were good
    % index = cellfun(@isempty, compResponses);
    c_resp{method+1} = compResponses;
    % resp(index(:, 1), :) = [];
    % strctCells_sigRamp(index(:, 1)) = [];
    % rL = cell2mat(resp(:, 2));
    % l_cells = strctCells_sigRamp(rL>350);
    % l_resp = resp(rL>350, :);
    
    % resp
    comparison(:, method+2) = compResponses(:, 2);
    
    
end
responses = compResponses;


% % checking if pruned cells were good
% strctCells_sigRamp = strctCells;
% index = cellfun(@isempty, compResponses);
% resp = compResponses;
% resp(index(:, 1), :) = [];
% strctCells_sigRamp(index(:, 1)) = [];
% rL = cell2mat(resp(:, 2));
% l_cells = strctCells_sigRamp(rL>350);
% l_resp = resp(rL>350, :);

% save([diskPath filesep taskPath filesep 'AllITCells_500stim_ReScreen.mat'], 'strctCells', 'psths', 'responses', '-v7.3')

%% for ease of viewing

pathStim = [diskPath filesep taskPath filesep 'forPaper' filesep 'AllCells' filesep 'SigandNonSig_Combined'];

imDir = dir(fullfile(pathStim));
imDir = imDir(~ismember({imDir.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));

% gets rid of weird copies
stimNames = struct2cell(imDir);
stimNames = stimNames(1, :)';
goodStim = ~startsWith(stimNames, '._', 'IgnoreCase', true);
imDir = imDir(goodStim);

[~, natIdx] = natsortfiles({imDir.name});


% imDir = imDir(natIdx);
%%

for i = 1:length(strctCells)
    bigNames{i, 1} = [num2str(strctCells(i).ChannelNumber) '_' num2str(strctCells(i).Name)];
end

%%
leftCells = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));

l_sc = bigNames(leftCells, :);
r_sc = bigNames(~leftCells, :);

[al, bl] = sort((l_sc(:, 1)));

[ar, br] = sort((r_sc(:, 1)));

l_comp = comparison(leftCells, :);
r_comp = comparison(~leftCells, :);
%%
sorted_comparison = [l_comp(bl, :); r_comp(br, :)];

%% histogram of response latencies (vis resp cells)

f = figure;
% histogram(cell2mat(sorted_comparison(:, 3)), 50:20:300, 'FaceColor', [0.6350 0.0780 0.1840]);
histogram(cell2mat(c_resp{1, 4}(:, 2)), 0:20:300, 'FaceColor', [0.6350 0.0780 0.1840]); % poisson latency

% stdErr = std(rL{4})/sqrt(length(rL{4}));
stdErr = std(rL{4});

% ylim([0 60])
title('Response latency of IT neurons')
xlabel('Response Latency (ms)');
ylabel('No of neurons');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
str = {'Mean = ', [num2str(mean(cell2mat(c_resp{1, 4}(:, 2))),' %.2f') ' +/- ' num2str(stdErr, ' %.2f')]};
text(200, 100, str, 'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold');
% xlim = [0 300];
% filename = [diskPath filesep taskPath filesep 'Hist_RespLat_Method1_2point5StdDev_80offset'];
filename = [diskPath filesep taskPath filesep 'Hist_SingleTrial_RespLat_Poisson']; % median 131.74ms
print(f, filename, '-dpng', '-r300')

%% histogram of ramp tuned neurons

pvals_all = cat(1, strctCells_sigRamp(:).pvalRamp);
sigInds = pvals_all < 0.01;
ramp_respLat = c_resp{1, 4}(sigInds, 2);

f = figure;
% histogram(cell2mat(sorted_comparison(:, 3)), 50:20:300, 'FaceColor', [0.6350 0.0780 0.1840]);
histogram(cell2mat(ramp_respLat), 50:20:300, 'FaceColor', [0.6350 0.0780 0.1840]);

stdErr = std(cell2mat(ramp_respLat))/sqrt(length(ramp_respLat));
% ylim([0 60])
title('Response latency of IT neurons')
xlabel('Response Latency (ms)');
ylabel('No of neurons');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
str = {'Mean = ', [num2str(mean(cell2mat(ramp_respLat)),' %.2f') ' +/- ' num2str(stdErr, ' %.2f')]};
text(200, 60, str, 'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold');

% filename = [diskPath filesep taskPath filesep 'Hist_RespLat_Method1_2point5StdDev_80offset'];
filename = [diskPath filesep taskPath filesep 'Hist_SingleTrial_RespLat_Poisson_SigRamp']; % median 130.585ms
% print(f, filename, '-dpng', '-r3000')

%% 
figure;

for i = 1:4
subplot(2, 2, i)
histogram(rL{i})
end

%% compute pval of ramp

% load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_alexnet_comparison_500Stimuli.mat']); % will create params 1x4 cell
% par = params;

% load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'scores_alexnet_comparison.mat']); % will make cell array with fc6 fc6 after relu fc7 fc7 after relu

layermat = 'fc6';
load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);


% params = score(:, 1:50);

netCells = [];
tic % takes 45 min to run with the p val distribution comp
for netwkNum = 1%:4
    
    
    %     params = par{netwkNum}(:, 1:50);
    
    pvalRamp = {};
    
    
    
    for i = 4%1:4%1:2 % testing the 3 different resp lat conputations
        
        strctCells_sigRamp = strctCells;
        
        %         if i == 1
        %             index = cellfun(@isempty, responses);
        %             resp = responses;
        %             resp(index(:, 1), :) = [];
        %             strctCells_sigRamp(index(:, 1)) = [];
        %         else
        %             index = cellfun(@isempty, compResponses);
        %             resp = compResponses;
        %             resp(index(:, 1), :) = [];
        %             strctCells_sigRamp(index(:, 1)) = [];
        %         end
        compResponses = c_resp{i};
        
        index = cellfun(@isempty, compResponses);
        resp = compResponses;
        resp(index(:, 1), :) = [];
        strctCells_sigRamp(index(:, 1)) = [];
        
        
        options.screenType = 'Object';
        
        for cellIndex = l(strctCells_sigRamp)
            %             n_reps = 100;
            options.ind_train = imageIDs; % use all objects to calculate STA
            %             p_dist = nan(1, n_reps);
            %             for dist = 1:n_reps
            dist = 1;
            [p, pRat] = Utilities.ObjectSpace.linearity_measure_STA(resp{cellIndex, 1}, params, options);
            %                 p_dist(dist) = p;
            %             end
            p_sum = p;
            %             % maybe median is a better way to summarize this? p_sum = median(p_dist)
            %             p_sum = sum(0.01 < p_dist)/length(p_dist);
            
            pvalRamp{cellIndex, i} = p_sum;
            strctCells_sigRamp(cellIndex).pvalRamp = p_sum;
%             disp(['Finished for cell ' num2str(cellIndex)])
        end
        
        pvals = cell2mat(pvalRamp(:, i));
        temp = [];
        for idx = 1:length(pvals)
            if pvals(idx) > 0.01
                temp = [temp idx];
            end
            
        end
        disp(length(strctCells_sigRamp) - length(temp));
        netCells(netwkNum, i) = length(strctCells_sigRamp) - length(temp);
        
    end
end
toc
%%
for i = 1:length(strctCells_sigRamp)
    
    if strctCells_sigRamp(i).pvalRamp > 0.01
        ind(i) = 1;
    end
    
end
ind = logical(ind);

strctCells = strctCells_sigRamp;
responses = resp;

strctCells(ind) = [];
responses(ind, :) = [];


%% carving up cells by category response

% get rid of non-responsive units
index = cellfun(@isempty, compResponses);
compResponses(index(:, 1), :) = [];
psths(index(:, 1), :) = [];
strctCells(index(:, 1)) = [];

cR = compResponses;
sig_cell = pvals < 0.014;
cR = compResponses(sig_cell, :);
strctCells_cR = strctCells(sig_cell);
strctCells_all = strctCells;
responses_cR = cR;

assert(isequal(responses_cR, cR));
psths_cR = psths(sig_cell);
% faceCells_ids = find(cell2mat(cR(:, 4)) == 1);
% textCells_ids = find(cell2mat(cR(:, 4)) == 2);
% vegCells_ids = find(cell2mat(cR(:, 4)) == 3);
% animCells_ids = find(cell2mat(cR(:, 4)) == 4);
% objCells_ids = find(cell2mat(cR(:, 4)) == 5);
%
%
% faceCells = cR(faceCells_ids, :);
% textCells = cR(textCells_ids, :);
% vegCells = cR(vegCells_ids, :);
% animCells = cR(animCells_ids, :);
% objCells = cR(objCells_ids, :);
%
% save([diskPath filesep 'Object_Screening' filesep 'ITResponses_SigRamp_SplitByCategory_500Stimuli'], 'faceCells', 'textCells', 'vegCells', 'animCells', 'objCells')
%%
for cat = 1:5
    
    ids = find(cell2mat(cR(:, 4)) == cat);
    strctCells = strctCells_cR(ids);
    responses = responses_cR(ids, :);
    psths = psths_cR(ids, :);
    
    switch cat
        case 1
            group = 'faceCells';
        case 2
            group = 'textCells';
        case 3
            group = 'vegCells';
        case 4
            group = 'animCells';
        case 5
            group = 'objCells';
    end
    
    save([diskPath filesep 'Object_Screening' filesep 'splitByCategory' filesep ['IT' group '_SigRamp_500Stimuli.mat']], 'strctCells', 'responses', 'psths')
    
end

%% scatterplot comparison of response latancies
% %
% index = cellfun(@isempty, compResponses);
% rL_resp = cell2mat(compResponses(:, 2));
% rL_resp(index(:, 1), 1) = 0; % only for plotting
% disp(sum(index(:, 1)))
%
%
% index = cellfun(@isempty, responses);
% rL_sel = cell2mat(compResponses(:, 2));
% rL_sel(index(:, 1), 1) = 0; % only for plotting
% disp(sum(index(:, 1)))
%
% f = figure; clf;
% hold on
% scatter(1:length(strctCells), rL_resp)
%
% scatter(1:length(strctCells), rL_sel)
%
%

