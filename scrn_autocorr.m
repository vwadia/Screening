
%% AUTOCORRELOGRAMS

% list of models one can try to have the curve fitting tookbox fit
% https://www.mathworks.com/help/curvefit/list-of-library-models-for-curve-and-surface-fitting.html

%% Recipes

% 1.
% - Take baseline raster (trials x time) and bin spikes (or not)
% - Then zscore each trial and compute xcorr(trial, trial) for all trials and average.
% - Fit per cell to see if there are cells with high RMSe (plot distribution and look)
% - Then average across cells in region and compute tau

% 2. Anna Levina's way - AC(t) = Cov[A(s), A(s+t)]/Var[A(s)]
% - Take baseline raster and convert to vector. Do this for all cells and then average those vectors
% resulting in a single 'multi-unit' train, average across all 8 channels to get 1 vector per region.
% - compute autocorr like this Cov[A(s), A(s+t)]/Var[A(s)] and loop over all t's (displacements between the two vectors - see image in 'Pictures' folder)
% Note: cov(A, B) % x = [A B], xc = x with mean of each column subtracted, cov = (xc' * xc)./n-1 where n = numel(A) = numel(B);
% Note: Cov[A(s), A(s+t)] --> cov of the overlap between the 2 where A(s+t) is a displaced version of A(s) by t timebins

%% CURRENT ISSUES 5/9/2022

% Tau values are not consistent across task or even different parts of the same task (Rec Enc vs Rec CR)
% MFC values are always extremely low (<100ms) 

% Need to fit per cell and plot histogram

% Combine different sessions (all cells) and run on large numbers of units 

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

binSpikes = 1; 
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


%% separate corr by area

rr = {};
ctr = 1;
for i = unique(areas') % vector has to be horizontal for i = vector to work in for loop
    
    if ~combineTraces
        rr{ctr} = R_t(areas == i, :);
    else
        rr{ctr} = R_t(ctr, :);
    end
    
    ctr = ctr+1;
    
end

[~, b] = sort(areas);
areaNames = areaNames(b);


%% plot autocorrs
arN = unique(areaNames, 'stable');
cols = Utilities.distinguishable_colors(length(arN));

xvals = length(binnedRas);

f1 = figure; clf;

hold on
title('Autocorrelations by area')
ct = 1;
for ar = 1:length(unique(areas))
   
    plot(mean(rr{ar}, 1, 'omitnan'), 'Color', cols(ct, :))
    ct = ct+1;
    
end
lgnd = legend(arN);

if binSpikes
    filename = [basePath filesep 'Autocorrelations_per_region_' task '_binWidth_' num2str(binWidth) '_method_' num2str(method) '_combineTraces_' num2str(combineTraces) '_trialWise_' num2str(trialWise)];
else
    filename = [basePath filesep 'Autocorrelations_per_region_' task '_method_' num2str(method) '_combineTraces_' num2str(combineTraces) '_trialWise_' num2str(trialWise)];
end
print(f1, filename, '-dpng', '-r0')

%% fit exp for timescales

% % % Juri's way
% modelfun   = @(b,x)(b(2)*(exp(-b(3)*x)+b(1)));
% % modelfun   = @(a,b,t,x)(a*(exp(-1/t*x)+b));
% beta0    = [0;0;0]; % initial guess
% opts     = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare'; % robust fitting
% opts.FunValCheck = 'off';
% 
% beta = {};
% tau_nl = [];
% for i = 1:length(unique(areas))
%     if exist('lags', 'var')
%         beta{i} = nlinfit(lags(lags>=0), mean(rr{i}, 1, 'omitnan'), modelfun, beta0, opts);
%     elseif exist('maxlag', 'var')
%         beta{i} = nlinfit([1:maxlag], mean(rr{i}, 1, 'omitnan'), modelfun, beta0, opts);
%     end
%     tau_nl(i) = 1/beta{i}(3);
% 
% 
% 
% end
% 
% figure;
% hold on
% title('tau\_nl')
% scatter(1:length(tau_nl), tau_nl, [], cols, 'filled');
% xticks([1:length(rr)])
% xticklabels([arN])

%% fit per cell

RMSe = {};
RMSe2 = {};
matlab_RMS = {};

tau = {};

for reg = 1:length(rr)
    
    region = rr{reg}; % all the cells in a region
    
    
    
    for i = 1:size(region, 1)
        
        if ~isnan(region(i, :))
            if exist('lags', 'var')
                ft{i} = fit(lags(lags>=0)', nanmean(region(i, :), 1)', 'exp1', 'Startpoint', [0, 0]);
                fitVec = ft{i}.a*exp(ft{i}.b*lags(lags>=0));
                
            elseif exist('maxlag', 'var')
                ft{i} = fit([1:maxlag]', nanmean(region(i, :), 1)', 'exp1', 'Startpoint', [0, 0]);
                fitVec = ft{i}.a*exp(ft{i}.b*[1:maxlag]);
                
            end
           
            tau{reg}(i) = (-1/ft{i}.b*1e3); % convert to ms (is this correct)?
            
            % RMSe
            realVec = nanmean(region(i, :), 1);
            RMSe{reg}(i) = sqrt(mean((realVec - fitVec).^2));
            RMSe2{reg}(i) = norm(realVec - fitVec);
        end
    end
    
end

if ~combineTraces
    m_tau = [];
    figure;
    hold on
    for r = 1:length(rr)
        histogram(RMSe{r});
        %         keyboard
        
        m_tau(r) = mean(tau{r});
        
    end
end
%% fit by region

tau_r = [];
RMSe_r = [];

for i = 1:length(rr)
    
    
    if exist('lags', 'var')
        ft_r{i} = fit(lags(lags>=0)', nanmean(rr{i}, 1)', 'exp1', 'Startpoint', [0, 0]);
        fitVec = ft_r{i}.a*exp(ft_r{i}.b*lags(lags>=0));
        
    elseif exist('maxlag', 'var')
        ft_r{i} = fit([1:maxlag]', nanmean(rr{i}, 1)', 'exp1', 'Startpoint', [0, 0]);
        fitVec = ft_r{i}.a*exp(ft_r{i}.b*[1:maxlag]);
        
    end
    
    tau_r(i) = (-1/ft_r{i}.b*1e3); % convert to ms (is this correct)?
    
    realVec = nanmean(rr{i}, 1);
    RMSe_r(i) = sqrt(mean((realVec - fitVec).^2));
    
    
end



f2 = figure;
hold on
title('tau per region')
scatter(1:length(tau_r), tau_r, [], cols, 'filled');
% ylim([0 500])
xticks([1:length(rr)])
xticklabels([arN])

if binSpikes
    filename = [basePath filesep 'Tau_per_region_' task '_binWidth_' num2str(binWidth) '_method_' num2str(method) '_combineTraces_' num2str(combineTraces) '_trialWise_' num2str(trialWise)];
else
    filename = [basePath filesep 'Tau_per_region_' task '_method_' num2str(method) '_combineTraces_' num2str(combineTraces) '_trialWise_' num2str(trialWise)];

end
print(f2, filename, '-dpng', '-r0')



%%  manual xcorr comp - understand how this works

% A = rand(100, 1);
% B = 2*A;
% res = zeros(1, size(A, 1) + size(B, 1) - 1);
% for i = 1:size(A, 1)+size(B, 1)-1
%     %     if i == size(A, 1)/2
%     %         keyboard
%     %     end
%     arg = (i - size(A, 1));
%     
%     if arg < 0
%         neg = 1;
%         limit = size(A, 1) + arg;
%     else
%         neg = 0;
%         limit = size(A, 1) - arg;
%     end
%     
%     for lim = 1:limit
%         if neg
%             res(i) = res(i) + A(lim) * B(lim - arg);
%             %             disp([lim, lim-arg])
%         else
%             res(i) = res(i) + A(lim + arg) * B(lim);
%             %             disp([lim+arg, lim])
%             
%         end
%         
%     end
%     
%     %     disp('--------');
%     
% end











