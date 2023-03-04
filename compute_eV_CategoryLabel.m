% Script to compute explained variance via category label 
% Do this for all responsive cells

% vwadia Nov2022

%% 

setDiskPaths 

taskPath = 'Object_Screening';
addpath(genpath('ObjectSpace'))

load([diskPath filesep 'Object_Screening' filesep 'AllMergedRespITCells_withPDist_Scrn_500Stim.mat']); RespCells = 1; layermat = 'fc6'; layerpy = 'IT_output'; stimDir = '500Stimuli'; 
% load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500Stim_Scrn_SigRamp.mat']); SigCells = 1; layermat = 'fc6'; layerpy = 'IT_output'; stimDir = '500Stimuli'; 

%% compute ev using category label and LOO

% make category labels for 500 ims
if strcmp(stimDir, '500Stimuli')
    imageIDs = [1:500]';
end
faceInds = 134:210;
objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
textInds = [264:282 290 400:408];
vegInds = [211:235 357:399];
animInds = [1:84 256:263];

labels = zeros(length(imageIDs), 1);
labels(ismember(imageIDs, faceInds)) = 1;
labels(ismember(imageIDs, textInds)) = 2;
labels(ismember(imageIDs, vegInds)) = 3;
labels(ismember(imageIDs, animInds)) = 4;
labels(ismember(imageIDs, objInds)) = 5;



pred_fr = nan(length(strctCells), 1);
ev_cat = nan(length(strctCells), 1);
ev_axis = nan(length(strctCells), 1);
% per cell 
for cellIndex = 1:length(strctCells)
    
    
    X = [labels ones(length(labels), 1)];
    
    Y = responses{cellIndex, 1};
    
    % use LOO
    for im = 1:length(Y)
        
        inds_train = setdiff(imageIDs, im);
        
        X_train = X(inds_train, :);
        
        Y_train = Y(inds_train);
        
        X_test = X(im, :);
        
        [b, ~, ~, ~, stats] = regress(Y_train, X_train);
        % stats(1) is the R^2 value
        
        pred_fr(im, 1) = X_test*b;
        
        
    end
    
    ev_cat(cellIndex) = Utilities.computeExplainedVariance(Y, pred_fr);
    
    
end

%% compute ev using axis model

load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);

for cellIndex = 1:length(strctCells)
    
    % collect responses
    resp = responses{cellIndex, 1};
    
    % calculate ev
    %             method = 'linear_regression';
    method = 'sta';
    [ev] = Utilities.ObjectSpace.STA_sub_cross_val(resp, params, method, 0.01);
    ev_axis(cellIndex, 1) = ev;
    
end

if strcmp(stimDir, '500Stimuli') && exist('SigCells')
    load([diskPath filesep taskPath filesep 'ExV_500Stim_1000Reps_SigRampCells']);
elseif strcmp(stimDir, '500Stimuli') && exist('RespCells')
    load([diskPath filesep taskPath filesep 'ExV_500Stim_1000Reps_RespCells']);    
end
expble_var = mean(explainableVar, 2);

%% compute ev using exemplar model - incomplete

% exemplar model - defined as a 3rd degree polynomial of the euclidean distance between exemplar and other face
% resp = c_0 + c_1*d + c_2*d^2 + c_3*d^3
%      = c_0 + c_1*(SUM_1toN(e_i - x_i))^0.5 + c_2*(SUM_1toN(e_i - x_i)) + c_3*(SUM_1toN(e_i - x_i))^1.5

load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);

for cellIndex = 1:llength(strctCells)
    
     % collect responses
    resp = responses{cellIndex, 1};
    
    fr = resp-mean(resp);
    
    % solve for the exemplar
    
    
    
end


%% noise ceiling 

A = cellfun(@(x) x', psths(:, 1), 'UniformOutput', false);

for cellIndex = 1:length(strctCells)
    ras = A{cellIndex, 1};
    if strcmp(strctCells(cellIndex).SessionID, 'P71CS_Fast')
        order = repelem(imageIDs, 6);
    elseif strcmp(strctCells(cellIndex).SessionID, 'P76CSRec_ReScreen')
        order3 = load([diskPath filesep 'Recall_Task' filesep 'P76CS' filesep 'ReScreenRecall_Session_1_20210917' filesep 'PsthandResponses']); % had a few trials cut off
        order = order3.order;
    else
        order = repelem(imageIDs, 4);
    end
    for im = 1:length(imageIDs)
        stimRas = ras(:, find(order == imageIDs(im)));
        A{cellIndex, im} = stimRas;
    end    
end
% ctr = 1;
for i=1:size(A,1)
    for j=1:size(A,2)
      hist=A{i,j}; 
      if strcmp(strctCells(i).SessionID, 'P71CS_Fast')
          stimDur = 167;
      else
          stimDur = 267;
      end
      if responses{i, 2}+stimDur > size(hist, 1)
          B{i,j}=mean(hist(responses{i, 2}:end,:),1);% Choose a time window you want to accumulate the spikes into a firing rate 
%           disp(i);
%           ctr = ctr + 1;
      else
          B{i,j}=mean(hist(responses{i, 2}:responses{i, 2}+stimDur,:),1);% Choose a time window you want to accumulate the spikes into a firing rate
      end
      fir(i,j)=mean(B{i,j});  %mean firing rate
    end
end
%compute noise ceiling
for i=1:size(B,1)
    for j=1:size(B,2)
        rast=B{i,j}';
        num(i,j)=length(rast);
        if length(rast)>1
            meu(i,j)=std(rast,0,1)/sqrt(length(rast));
        end
    end
    ev00(i)=1-sum(meu(i,num(i,:)>1).^2)/sum((fir(i,num(i,:)>1)-mean(fir(i,num(i,:)>1))).^2);
end
evn=ev00;

%%
highnc = (evn > 0.1)';

% tokeep = (highnc + highev)  == 2;
tokeep = highnc; 

expble_var = expble_var(tokeep);
ev_cat = ev_cat(tokeep);
ev_axis = ev_axis(tokeep);

keptCells = strctCells(tokeep);
%% figure 



% h1 = ev_cat./expble_var;
% h2 = ev_axis./expble_var;
% 
% h1(h1 > 1) = 1;
% h2(h2 > 1) = 1;
% % histogram of ev/exv
% f = figure; 
% hold on; 
% histogram(h1, -0.2:0.1:1, 'FaceColor', [0.6350 0.0780 0.1840], 'EdgeColor', [0 0 0]); 
% histogram(h2, -0.2:0.1:1, 'FaceColor', [0.4940 0.1840 0.5560], 'EdgeColor', [0 0 0]);
% 
% title({'Explained variance of IT neurons', 'Category label vs Axis model'}, 'FontSize', 16, 'FontWeight', 'bold')
% xlabel('Explained variance')
% ylabel('Number of cells')
% lgnd = legend({'Category label', 'Axis model'}, 'FontWeight', 'bold');
% 
% set(gca, 'FontSize', 14, 'FontWeight', 'bold');
% 
% filename = [diskPath filesep taskPath filesep 'Hist_eV_' method 'vsCatLab _forPaper'];
% 
% if exist('RespCells', 'var') && RespCells == 1
% 
%     filename = [filename '_RespCells'];
% elseif exist('SigCells', 'var') && SigCells == 1
% 
%     filename = [filename '_SigRampCells'];
% end
% 
% print(f, filename, '-dpng', '-r0')
% close all


%% cdf 

f = figure;
hold on
h1 = ev_cat./expble_var;
h2 = ev_axis./expble_var;


% pvals_ramp = pvals_ramp(tokeep);
h1(h1 >= 1) = 1;
h2(h2 >= 1) = 1;

h1(h1 <= -0.2) = [];
h2(h2 <= -0.2) = [];

[rej_null, p] = kstest2(h1, h2);

cd1 = cdfplot(h1);
cd2 = cdfplot(h2);
cd1.LineWidth = 2;
cd1.Color = [0.6350 0.0780 0.1840];
cd2.LineWidth = 2;
cd2.Color = [0.4940 0.1840 0.5560];

xlim([-0.2 1])

x_lim = xlim;
y_lim = ylim;

title({'Explained variance of visually responsive IT neurons', 'Category label vs Axis model'}, 'FontSize', 16, 'FontWeight', 'bold')
xlabel('x = Explained variance')
ylabel('F(x)')
lgnd = legend({'Category label', 'Axis model'}, 'FontWeight', 'bold');
lgnd.Position = [0.692,0.631,0.1848,0.0631];
text(x_lim(1)*0.80, y_lim(2)*0.88, ['p = ' num2str(p)], 'FontSize', 14,'FontWeight', 'bold');

set(gca, 'FontSize', 14, 'FontWeight', 'bold');

filename = [diskPath filesep taskPath filesep 'CDF_eV_' method 'vsCatLab_forPaper'];

if exist('RespCells', 'var') && RespCells == 1

    filename = [filename '_RespCells'];
elseif exist('SigCells', 'var') && SigCells == 1

    filename = [filename '_SigRampCells'];
end

print(f, filename, '-dpng', '-r0')
