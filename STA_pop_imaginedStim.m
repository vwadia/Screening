%% STA imagination work
%% load in data

setDiskPaths
% task = 'Object_Screening';
task = 'Recall_Task';
taskPath = [diskPath filesep task filesep 'forPaper'];
if ~exist(taskPath)
    mkdir([taskPath]);
end


addpath(genpath('osortTextUI'));
addpath(genpath('ObjectSpace'));
addpath(genpath('synthetic_face_generator'));
% how is binned average currently computed?
    % Are the computed 'scresponses' to each image normalized? - if not need
    % to do this?
% What are the xy positions in PC space of the recalled stim for each cell?
% - need to save these and their corscresponding FRs
% load([diskPath filesep 'Recall_Task' filesep 'AllCells_500stim_Im']);
% load([diskPath filesep 'Recall_Task' filesep 'Allscresponses_500stim_Im']);

if strcmp(task, 'Recall_Task')
    load([diskPath filesep 'Recall_Task' filesep 'AllITCells_500stim_Im_SigRamp']);
    load([diskPath filesep 'Recall_Task' filesep 'AllITResponses_500stim_Im_SigRamp']);
    lbl = 'Imagination';

elseif strcmp(task, 'Object_Screening')
    load([diskPath filesep 'Object_Screening' filesep 'MergedITCells_500stim_Scrn_SigRamp']); % loads cells with sig ramps
    lbl = 'Viewing';
    
end
load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']); % will create params = 500x50

separateSides = 0;
% % per patientc
% pt_ID = 'P79CS';
% % pt_ID = 'P76CS';
% 
% for i = l(strctCells)
%     ids_toKeep(i, 1) = strcmp(strctCells(i).PatientID, pt_ID);
%     
% end
% strctCells = strctCells(ids_toKeep);
% strctResp = strctResp(ids_toKeep);
% responses = responses(ids_toKeep, :);
% psths = psths(ids_toKeep, :);

%% separating hemispheres
if separateSides
    left = 0;
    right = 1;
    strctCELL = struct2cell(strctCells');
    strctCELL = strctCELL';
    
    if left
        ITSide = cellfun(@(x) strcmp(x, 'LFFA'), strctCELL(:, 4));
        
    elseif right
        ITSide = cellfun(@(x) strcmp(x, 'RFFA'), strctCELL(:, 4));
        
    end
    
    strctCells = strctCells(ITSide);
    psths = psths(ITSide, :);
    responses = responses(ITSide, :);
    if exist('strctResp', 'var')
        strctResp = strctResp(ITSide);
    end

end
%% population STA plot in imagination
xvals = [];
yvals = [];
if strcmp(task, 'Recall_Task')
    % concatenate frs to all ims across cells (mean subbed) and all positions in both x and y - correlate them
    options.nbins = 4;
    
    Im_resp = [];
    for cellIndex = l(strctCells)
        xvals = [xvals strctCells(cellIndex).Im_xvals];
        yvals = [yvals; strctCells(cellIndex).Im_yvals];
        
        cellresp = strctResp(cellIndex).CRResp;
        %     Im_resp = [Im_resp; (cellresp - mean(cellresp))./std(cellresp)];
        Im_resp = [Im_resp; (cellresp - mean(cellresp))];
        
    end
    xvals = xvals';
    touse_resp = Im_resp;
    disp('done')
elseif strcmp(task, 'Object_Screening')
    options.screenType = 'Object';
    options.ind_train = [1:500]';
    s_resp = [];
    for cellIndex = l(strctCells)
        
        cellresp = responses{cellIndex, 1};

        vals = Utilities.ObjectSpace.grab_AxisProj_Values(responses{cellIndex, 1}, params, options);
        xvals = [xvals; vals(:, 1)];
        yvals = [yvals; vals(:, 2)];
        
        
        s_resp = [s_resp; (cellresp - mean(cellresp))];
        
    end
   touse_resp = s_resp; 
end
% shuffle control
n_reps = 1000;
cc_rand = zeros(n_reps, 1);
cc_ortho_rand = zeros(n_reps, 1);
for i = 1:n_reps
    
    cc_rand(i) = corr(Utilities.Shuffle(xvals), touse_resp);
    cc_ortho_rand(i) = corr(Utilities.Shuffle(yvals), touse_resp);
end
cc = corr(xvals,touse_resp); % 
p = sum(cc_rand > cc)/length(cc_rand);
cc_ortho = corr(yvals, touse_resp);
p_ortho = sum(cc_ortho_rand > cc_ortho)/length(cc_ortho_rand);


%% make figure
f = figure;
sgtitle({'Average projections onto individual axes', 'Imagination - all cells'})
% FR along "STA"
h1 = subplot(3, 3, [2 3])
nonlin = Utilities.ObjectSpace.compute_binned_average(xvals, Im_resp, options.nbins, 10);
errorbar(nonlin.x, nonlin.y, nonlin.e, 'k');
yl_sta = ylim;
xl_sta = xlim;
% ylim([-yl_sta(2) yl_sta(2)])

% scatter
dot_color = zeros(size(Im_resp,1),3); 
dot_color(:,1) = ((Im_resp-min(Im_resp))/(max(Im_resp)-min(Im_resp)));
dot_color(:,3) = 1- dot_color(:,1);

[~, reord] = sort(Im_resp);
xvals_s = xvals(reord);
yvals_s = yvals(reord);
c = dot_color(reord, :);
h2 = subplot(3, 3, [5 6 8 9])
scatter(xvals_s, yvals_s, 30, c, 'filled');
box on

% FR along ortho
h3 = subplot(3, 3, [4 7])
nonlin = Utilities.ObjectSpace.compute_binned_average(yvals, Im_resp, options.nbins, 10);
herrorbar(nonlin.y, nonlin.x, nonlin.e, 'k');
xlim([-yl_sta(2) yl_sta(2)]);
ylim([xl_sta(1) xl_sta(2)]);
% shuffle dist
h4 = subplot(3, 3, 1)
h = histogram(cc_rand,-0.5:0.01:1);
h.FaceColor = [1 1 1];
hold on;
plot([cc cc],[0 150], 'LineWidth', 2, 'Color', 'r');
text(.2, 200, num2str(p,'p = %.3f'))

linkaxes([h1, h2], 'x');
linkaxes([h2, h3], 'y');

filename = [taskPath filesep 'RampSummaryAxisProj_AllCells_Im'];
% print(f, filename, '-dpng', '-r0');

%% Shuffle distribution along each axis - imagination

% ortho
f = figure; hold on; 
h = histogram(cc_ortho_rand, -1:0.1:1); 
h.FaceColor = [1 1 1]; 
plot([cc_ortho cc_ortho], [0 500], 'LineWidth', 2, 'Color', 'r'); 
text(.3, 1000, num2str(p_ortho, 'p = %.3f'))

if strcmp(task, 'Recall_Task')
    title({'Correlation of imagined responses and projection value', 'Orthogonal axis'})
elseif  strcmp(task, 'Object_Screening')
    title({'Correlation of viewed responses and projection value', 'Orthogonal axis'})
end
ylabel('Bin count');
xlabel('Correlation value'); 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
filename = [taskPath filesep 'RampSummaryAxisProj_AllCells_OrthoAxCorr']
filename = [filename '_' lbl];
% print(f, filename, '-dpng', '-r0');



% sta
f = figure; hold on;
h = histogram(cc_rand,-1:0.01:1);
h.FaceColor = [1 1 1];
hold on;
if strcmp(task, 'Recall_Task')
    plot([cc cc],[0 50], 'LineWidth', 2, 'Color', 'r');
    text(.2, 100, num2str(p,'p = %.3f'))
    title({'Correlation of imagined responses and projection value', 'Preferred axis'})
    
elseif strcmp(task, 'Object_Screening')
    plot([cc cc],[0 500], 'LineWidth', 2, 'Color', 'r');
    text(.2, 1000, num2str(p,'p = %.3f'))
    title({'Correlation of viewed responses and projection value', 'Preferred axis'})
    
end
ylabel('Bin count');
xlabel('Correlation value'); 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

filename = [taskPath filesep 'RampSummaryAxisProj_AllCells_STAAxCorr'];
filename = [filename '_' lbl];

% print(f, filename, '-dpng', '-r0');
%% comparing real and imagined scresponses along the ramp and ortho axes per cell

% axis = 'STA';
axis = 'ortho';

if strcmp(axis, 'STA')
    sf_r = zeros(length(strctCells), 8);    
    max_sf_r = zeros(length(strctCells), 1);
elseif strcmp(axis, 'ortho')
    sf_r_o = zeros(length(strctCells), 8);
end

for cellIndex = l(strctCells)
     if strcmp(axis, 'STA')
        crds = strctCells(cellIndex).Im_xvals;
        
        [~, idx] = sort(crds);
        if strctCells(cellIndex).SessionID == 2
            scresp = strctResp(cellIndex).EncResp';
        else
            scresp = strctResp(cellIndex).ScrnResp';
        end
        scresp = scresp(idx);
        scresp = scresp - mean(scresp);
        max_sf_r(cellIndex) = max(scresp);
        if cellIndex < 10
            sf_r(cellIndex, 3:end) = scresp./max_sf_r(cellIndex);           
        else
            sf_r(cellIndex, :) = scresp./max_sf_r(cellIndex);
        end
        
    elseif strcmp(axis, 'ortho')
        assert(length(unique(max_sf_r)) > 1, 'Need to compute max for normalization');
        crds = strctCells(cellIndex).Im_yvals';
        [~, idx] = sort(crds);
        if strctCells(cellIndex).SessionID == 2
            scresp = strctResp(cellIndex).EncResp';
        else
            scresp = strctResp(cellIndex).ScrnResp';
        end
        scresp = scresp(idx);
        scresp = scresp - mean(scresp);

        if cellIndex < 10
            sf_r_o(cellIndex, 3:end) = scresp./max_sf_r(cellIndex);           
        else
            sf_r_o(cellIndex, :) = scresp./max_sf_r(cellIndex);
        end
        
    end
    
end

%% compute correlations

cc_task_rand = zeros(n_reps, length(strctCells));
cc_task_ortho_rand = zeros(n_reps, length(strctCells));
cc_task = zeros(1, length(strctCells));
cc_task_ortho = zeros(1, length(strctCells));

for cellIndex = l(strctCells)
    crds_x = strctCells(cellIndex).Im_xvals;
    crds_y = strctCells(cellIndex).Im_yvals';
    [~, idx] = sort(crds_x);
    [~, idy] = sort(crds_y);
    
    fr = strctResp(cellIndex).CRResp - mean(strctResp(cellIndex).CRResp);
    fr = fr./std(fr);
    if strctCells(cellIndex).SessionID == 2
        sfr = strctResp(cellIndex).EncResp - mean(strctResp(cellIndex).EncResp);
    else
        sfr = strctResp(cellIndex).ScrnResp - mean(strctResp(cellIndex).ScrnResp);
    end
    sfr = sfr./std(sfr);
    
    fr = fr(idx);
    sfr = sfr(idx);
    
    fro = fr(idy);
    sfro = sfr(idy);

%     sfr = sf_r(cellIndex, :)';
%     fr = f_r(cellIndex, :)';
%     sfro = sf_r_o(cellIndex, :)';
%     fro = f_r_o(cellIndex, :)';
    
    % these should be the same as the values are simply shuffled
    cc_task(cellIndex) = corr(sfr, fr);
    cc_task_ortho(cellIndex) = corr(sfro, fro);
    
    for i = 1:n_reps
        cc_task_rand(i, cellIndex) = corr(Utilities.Shuffle(sfr), Utilities.Shuffle(fr));
        cc_task_ortho_rand(i, cellIndex) = corr(Utilities.Shuffle(sfro), Utilities.Shuffle(fro));
    end
    
    p(cellIndex) = sum(cc_task_rand(:, cellIndex) > cc_task(cellIndex))/length(cc_task_rand(:, cellIndex));
    p_ortho(cellIndex) = sum(cc_task_ortho_rand(:, cellIndex) > cc_task_ortho(cellIndex))/length(cc_task_ortho_rand(:, cellIndex));

end


%% plot correlations between viewed and Im responses per cell


% f = figure;
% set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
% 
% for cellIndex = l(strctCells)
%     s = subplot(3, 4, cellIndex);
%     h = histogram(cc_task_rand(:, cellIndex), -1:0.1:1);
%     h.FaceColor = [1 1 1];
%     hold on;
%     plot([cc_task(cellIndex) cc_task(cellIndex)],[0 50], 'LineWidth', 2, 'Color', 'r');
%     text(.2, 100, num2str(p(cellIndex),'p = %.3f'))
%     
% %     s = subplot(3, 4, cellIndex);
% %     h = histogram(cc_task_ortho_rand(:, cellIndex), -1:0.1:1);
% %     h.FaceColor = [1 1 1];
% %     hold on;
% %     plot([cc_task_ortho(cellIndex) cc_task_ortho(cellIndex)],[0 50], 'LineWidth', 2, 'Color', 'r');
% %     text(.2, 100, num2str(p_ortho(cellIndex),'p = %.3f'))
%     
% end
% sgtitle('Viewed-Imagined correlations per cell')
% filename = [taskPath filesep 'ViewedIm_correlations_perCell'];
% print(f, filename, '-dpng', '-r0');

% % Average across cells
% f = figure;
% set(gcf,'Position',get(0,'Screensize')) % display fullsize on other screen
% 
% p_test = sum(mean(cc_task_rand, 2) > mean(cc_task))/length(mean(cc_task_rand, 2));
% 
% h = histogram(mean(cc_task_rand, 2), -1:0.01:1);
% h.FaceColor = [1 1 1];
% hold on;
% plot([mean(cc_task) mean(cc_task)],[0 100], 'LineWidth', 2, 'Color', 'r');
% text(.6, 100, num2str(p_test,'p = %.3f'))
% 
% sgtitle('Viewed-Imagined correlation averaged across cells')
% filename = [taskPath filesep 'ViewedIm_correlation_avg'];
% print(f, filename, '-dpng', '-r0');



