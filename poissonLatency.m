% Script messing with Poisson GLM used by Brooks
% 
% vwadia Feb2023

%% 

setDiskPaths

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

timelimits = [-0.17 0.53]; %sec
%% separate cells by session (start with P73)

% sessCells = cellfun(@(x) strcmp(x, 'P73CS_ParamObj'), strctCELL(:, 8), 'UniformOutput', false);
% 
% sessCells = cell2mat(sessCells);
% 
% strctCells = strctCells(sessCells);
% psths = psths(sessCells, :);
% responses = responses(sessCells, :);

%% for response latency testing - cells with weird psths that this alg should get 

weirdCells = [2150 4516 949 1056 1383 661 713 1905 204 1755 3178 625 651 1036 1981 ];
allCells = cell2mat(strctCELL(:, 1));

weirdIds = find(allCells == weirdCells);
weirdIds = mod(weirdIds, length(strctCells));

weirdIds(weirdIds == 0) = length(strctCells); 

strctCells = strctCells(weirdIds); 
psths = psths(weirdIds, :);
responses = responses(weirdIds, :);

%%
% Open questions:
%     What is anchor? - used to compute cISI which is fed into poisscdf
%     What is the best time window to use for most accurate burst detection?

% Notes:
%     Significance and min spike # per burst is set inside p_burst
%     Function doesn't handle negative numbers

for cellIndex = 1:length(strctCells)
exCell = psths{cellIndex, 1}; % nice animal cell psth

% [BOB, EOB, SOB]=p_burst(InTrain, StartT, StopT,doDisplay,varargin)
% [m1, p1] = max(sum(exCell, 2));
BOB = zeros(1, size(exCell, 1));
EOB = zeros(1, size(exCell, 1));
SOB = zeros(1, size(exCell, 1));
stamps = zeros(size(exCell, 1), 3);
onTimes = [];

% avgSpikRate = mean(mean(exCell(:, -timelimits(1)*1e3:-timelimits(1)*1e3+50))); % baseline FR of cell
for it = 1:size(exCell, 1)
    
   
    % spike timestamps - note that function gets rid of negative numbers
    times = find(exCell(it, :) == 1);
    
    % how does one assess significance? - it's inside the function
    % how does one include a given cell's baseline? - start and stop time are used to compute average FR
    % note start time can't be 0
    startT = 170;
    endT = 170+267;
%     [b, e, s] = Utilities.p_burst(times, startT, endT, 0); 

%     can manually input an average spike rate - per trial works best so far
    avgSpikRate = sum(times > -timelimits(1)*1e3 & times < -timelimits(1)*1e3+50)/50; % baseline FR per trial
%     avgSpikRate = sum(times > -timelimits(1)*1e3)/-timelimits(1)*1e3; % baseline FR per trial
    [b, e, s] = Utilities.p_burst(times, startT, endT, 0, avgSpikRate); 
    
    if ~isempty(b)   
%         train = times(times > -timelimits(1)*1e3);
        train = times;%(times > -timelimits(1)*1e3);
        
        % in cases with multiple burtst take one with max surprise
        if length(s) > 1
            [~, pos] = max(s);
            stamps(it, 1) = train(b(pos));
            stamps(it, 2) = train(e(pos));
            stamps(it, 3) = s(pos);
        else % else take the first 
            stamps(it, 1) = train(b);
            stamps(it, 2) = train(e);
            stamps(it, 3) = s;
        end
    end
end

numTrials(cellIndex, 1) = sum(stamps(:, 1) ~= 0);
onTimes = stamps(find(stamps(:, 1) ~= 0), 1);
adjOnTimes = onTimes(onTimes > 170);
adjOnTimes = adjOnTimes - 170;
adj{cellIndex, 1} = adjOnTimes;
% train(BOB)
% train(EOB)
end

respLat(:, 1) = cellfun(@(x) mean(x), adj, 'UniformOutput', false);
respLat(:, 2) = cellfun(@(x) std(x), adj, 'UniformOutput', false);
