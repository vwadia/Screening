
%% Computing autocorrelation per cell
% 
% Keeping only those cells that show exponential fitting
%     Compute RMSe for each cell, view histogram, and create a cutoff per region
%     Do not use first time step - use step 2-4 with sharpest dropoff of autocorr

%%  

dbstop if error
[~, host] = system('hostname');
if strcmp(host(1:end-1), 'DWA644201')
    atCedars = 1;
    diskPath = 'G:\SUAnalysis';
elseif strcmp(host(1:end-1), 'DESKTOP-LJHLIED')
    atCedars = 0;
    diskPath = 'G:\SUAnalysis';
elseif strcmp(host(1:end-1), 'Varuns-iMac-2.local')
    atCedars = 0;
    diskPath = '/Volumes/T7/SUAnalysis';
end
taskPath = 'Object_Screening';

%% load in data

% screeningScript % set paths appropriately
timelimits = screeningData.timelimits; % this depends on the data loaded in
psth = screeningData.psth;

% % All cells in object Screen
% load([diskPath filesep taskPath filesep 'AllCells_500stim_Scrn']);
% % IT cells only
% % load([diskPath filesep taskPath filesep 'AllITCells_500stim_Scrn']);
% timelimits = [-0.13 0.53]; % this depends on the data loaded in
% psth = psths;
% task = 'Scrn';
% basePath = [diskPath filesep taskPath];

% RecallScript % set paths appropriately
% timelimits = [-RecallData.offsetEnc(1)*1e-3 RecallData.offsetEnc(2)*1e-3]; % this depends on the data loaded in
% 
% psth = RecallData.CRTimeCourse; 
% task = 'RecCR';


% psth = RecallData.EncodingTimeCourse; 
% task = 'RecEnc';

%% preprocess data

method = 1; % 1 - xcorr, 2 - Levina method
combineTraces = 0; % whether to create multiunit traces or keep raw signal
trialWise = 1; % compute trialwise autocorr and average VS creating a single vector per cell and computing autocorr of that

% Juri did 1 0 1 with LfP, Levina does 2 1 0

binSpikes = 0; 
binWidth = 50; % ms
dataMatrix = [];
tic
for cellIndex = 1:length(strctCells)
    
    pt = psth{cellIndex, 3} > -0.001 & psth{cellIndex, 3} < 0.001;
    startPoint = find(pt == 1);
        
    b_raster = psth{cellIndex, 1}(:, startPoint+(timelimits*1e3):startPoint); 
    
    if binSpikes
        binnedRas = [];       
        for bin = 1:binWidth:size(b_raster, 2) % so N_j stays in range
            if bin+binWidth-1 > size(b_raster, 2)
                binnedRas = [binnedRas mean(b_raster(:, bin:end), 2)];  
            else
                binnedRas = [binnedRas mean(b_raster(:, bin:bin+binWidth-1), 2)];         
            end
        end
    else
        binnedRas = b_raster;
    end
    
    if trialWise
        dataMatrix = cat(3, dataMatrix, binnedRas);
    else
        dataMatrix = [dataMatrix; binnedRas(:)'];
    end
  
    % brain area info
    areas(cellIndex, 1) = strctCells(cellIndex).brainAreaIndex;
    areaNames{cellIndex,1} = strctCells(cellIndex).brainArea;
    
end

% making MU traces per region
if combineTraces

    MU_trace = [];
    for i = unique(areas')
        if trialWise
            MU_trace = cat(3, MU_trace, mean(dataMatrix(:, :, areas == i), 3));
        else
            MU_trace = [MU_trace; mean(dataMatrix(areas == i, :), 1)]; % verage across all cells in the region
        end
    end
    
    dataMatrix = MU_trace;
    
end
toc

%% Compute autocorr

AC = [];
R_t = [];
clearvars lags maxlag;

if trialWise
    
    for cellIndex = 1:size(dataMatrix, 3)     
        
        R = [];
        
        ctr = 1;
        AC = [];
        
        for tr = 1:size(dataMatrix, 1)
            if sum(dataMatrix(tr, :, cellIndex)) ~= 0
                bR_trial = dataMatrix(tr, :, cellIndex);
                
                if method == 1
                    
                    bR_trial = zscore(bR_trial, [], 2); % zscoring each row

                    [R, lags] = xcorr(bR_trial, bR_trial, 'coeff'); % coeff makes sure autocorr at lag 0 == 1
                    AC(ctr, :) = R(lags > 0);
                    
                elseif method == 2
                    if binSpikes
                        maxlag = ceil(size(bR_trial, 2)/2);
                    else
                        maxlag = (-timelimits(1)*1e3)/2; 
                    end   
                    for t = 1:maxlag
                        cv = cov(bR_trial(t:end), bR_trial(1:end-t+1));
                        AC(ctr, t) = cv(1, 2)/var(bR_trial);
                    end
                    
                end
                ctr = ctr + 1;
            end
        end
        
        if isempty(AC)
            if exist('maxlag', 'var')
                AC = nan(1, maxlag);
            else
                AC = nan(1, size(bR_trial, 2)-1);
            end
        end
        
        R_t(cellIndex, :) = nanmean(AC, 1);
        
    end
    
else
    
    for cellIndex = 1:size(dataMatrix, 1)  
        
        bR_trial = dataMatrix(cellIndex, :);

        if method == 1
            bR_trial = zscore(bR_trial, [], 2);
            
            R = [];
            [R, lags] = xcorr(bR_trial, bR_trial, 'coeff'); % coeff makes sure autocorr at lag 0 == 1
            AC(cellIndex, :) = R(lags >= 0);
            
        elseif method == 2

            if binSpikes
                maxlag = ceil(size(bR_trial, 2)/2);
            else
                maxlag = (-timelimits(1)*1e3)/2;
            end

            for t = 1:maxlag
                cv = cov(bR_trial(t:end), bR_trial(1:end-t+1));
                AC(cellIndex, t) = cv(1, 2)/var(bR_trial);
            end
            
        end
        
    end  
    
    R_t = AC;

end

%% Why is this not working suddenly? won't work if you use x_r

cellIdx = 4;
[~, t] = min(diff(R_t(cellIdx, 1:20)));
% t = 13;
x = lags(lags>t & lags<151)';
x_r = linspace(1, length(x), length(x))'; % for plotting
y = nanmean(R_t(cellIdx, t+1:150), 1)';
assert(isequal(length(x), length(y)));


% % Juri's way
modelfun   = @(b,x)(b(2)*(exp(-b(3)*x)+b(1)));
beta0    = [-1;1;1]; % initial guess - [0;0;0] doesn't work
opts     = statset('nlinfit');
opts.RobustWgtFun = 'bisquare'; % robust fitting
opts.FunValCheck = 'off';


ft_nl = nlinfit(x, y, modelfun, beta0, opts);

fitVec = ft_nl(2)*(exp(-ft_nl(3)*x)+ft_nl(1));

realVec = y;
RMSe_r = sqrt(mean((realVec - fitVec).^2));

% plot and look
figure; hold on; scatter(x, y); plot(fitVec)


%%

f = figure;
hold on
for cellIndex = l(strctCells)
    
plot(R_t(cellIndex, :), '-o')
keyboard

end
