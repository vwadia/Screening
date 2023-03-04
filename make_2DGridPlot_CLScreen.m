
% Script to plot 2D grid of synthetic pictures for Synth screening
% for reference see figure 5 of Stevens paper

% my images are labeled cellidxyz_%03d_Grid_(xcrd,ycrd)
% so read in name, then x value = ( to , and y value = , to )
% save the vals and the corresponding fr
% make colors 
% draw 2D scatter plot with appropriate colors (exatly like STA plot)
% also make 2D STA plot with 4096D vector and see what is looks like 

% vwadiaJan2023


setDiskPaths

taskPath = 'Object_Screening';
patID = 'P82CS';
m_sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115'];
a_sessPath = [diskPath filesep taskPath filesep 'P82CS' filesep 'ClosedLoopReScreen_Session_1_20230115'];


stimPath = [a_sessPath filesep 'genStimOnly'];
    
load([a_sessPath filesep 'SynthPsthandResponses'])
a_strctCells = load([a_sessPath filesep 'strctCells']);

load([a_sessPath filesep 'PsthandResponses'])


m_strctCells = load([m_sessPath filesep 'strctCells']);

%% 

if strcmp(patID, 'P82CS')
    corrThresh = 0.9;
    morn_cells = [904];
    aft_matches = [2360];
    aft_idx = [1];
    
    n_steps = 5;
    n_steps_ortho = 5; % set per session
    scale = 2;
    stepRangeOrtho = [-n_steps_ortho:1:n_steps_ortho]*scale;
    stepRange = [-n_steps:1:n_steps]*scale;
    nbin = 10;
    
    greyImsOnly = 1; % only black and white images
end

fullImDir = Utilities.readInFiles(stimPath);
fullImDirCell = struct2cell(fullImDir)';

for cellIndex = 1:length(aft_matches)
    
%     mornCellIdx = structfind(m_strctCells.strctCells, 'Name', morn_cells(cellIndex));

    mornCell = morn_cells(cellIndex);
    
    % images corresponding to that cell
    cellIms = cellfun(@(x) strcmp(x(7:7+numel(num2str(mornCell))-1), num2str(mornCell)), fullImDirCell(:, 1), 'UniformOutput', false);
    
    idx = structfind(a_strctCells.strctCells, 'Name', aft_matches(cellIndex));
    
    if ~isempty(synthResponses{idx, 1})
        
        resp = synthResponses{idx, 1}(cell2mat(cellIms), 1); % responses only to desired images
          
        imNames = fullImDirCell(cell2mat(cellIms), 1); 
        
        % extract grid ims only
        % note I'm searching imNames because those are already images for
        % just that cell and going from dashpos+1 to dashpos+4 because I'm
        % not recomputing dashpos for every name
        dashpos = strfind(imNames{1}, '_');       
        cellGridIms = cellfun(@(x) strcmp(x(dashpos(1)+1:dashpos(1)+4), 'Grid'), imNames, 'UniformOutput', false); 
        
        imNames = imNames(cell2mat(cellGridIms));
         % only responses to grid images (for now)
        resp = resp(cell2mat(cellGridIms)); 
        
        if strcmp(patID, 'P82CS')
            if greyImsOnly
                % all the b&W images were first 
                imNames = imNames(1:120);
                resp = resp(1:120);            
            else
                imNames = imNames(121:end);
                resp = resp(121:end);
            end
        end
        
        for im = 1:length(imNames)
            
            brapos = strfind(imNames{im}, '(');
            ketpos = strfind(imNames{im}, ')');
            compos = strfind(imNames{im}, ',');
            
            x(im) = str2num(imNames{im}(brapos+1:compos-1));
            y(im) = str2num(imNames{im}(compos+1:ketpos-1));
            
        end
        
        
        
        [sorted_resp, reorder_ind] = sort(resp);
            
        f = figure;
        hold on
        
        if greyImsOnly
            sgtitle({'2D spread of synthetic images ', [patID ' B&W Images']}) 
        else
            sgtitle({'2D spread of synthetic images ', [patID ' Colored Images' ]}) 
        end
        % make plot 
        h2 = subplot(3, 3, [5 6 8 9]);
 
        minResp = min(resp);
        maxResp = max(resp);
        dot_color = zeros(size(resp,1),3);
        dot_color(:,1) = ((resp-minResp)/(maxResp-minResp));
        dot_color(:,3) = 1- dot_color(:,1);
        dot_size = 20;
        
        x = x(reorder_ind);
        y = y(reorder_ind);
        
        dot_color = dot_color(reorder_ind, :);
        
        scatter(x, y, dot_size, dot_color, 'filled')
        
        % binned FR along pref ax for all ortho values
        h1 = subplot(3, 3, [2 3]);
        hold on
        for nc = 1:length(stepRangeOrtho)
            imseq = find(y == stepRangeOrtho(nc));
            im_x = x(imseq);
            fr = sorted_resp(imseq);
            nonlin = compute_binned_average(im_x, fr, nbin, 1); 
            hline = plot(nonlin.x, nonlin.y);
            hline.Color = [0.5 0.5 0.5 0.5];  % alpha=0.5
        end
        
        % summary in black
        nonlin = compute_binned_average(x, sorted_resp, nbin, 1); 
        errorbar(nonlin.x, nonlin.y, nonlin.e, 'k');

        
        % same for ortho ax
        h3 = subplot(3, 3, [4 7]);
        hold on
        for nc = 1:length(stepRange)
            imseq = find(x == stepRange(nc));
            im_y = y(imseq);
            fr = sorted_resp(imseq);
            nonlin = compute_binned_average(im_y, fr, nbin, 1);
            hline = plot(nonlin.y - mean(nonlin.y), nonlin.x);
            hline.Color = [0.5 0.5 0.5 0.5];  % alpha=0.5
        end
        
        % summary in black
        nonlin = compute_binned_average(y, sorted_resp, nbin, 1); 
        herrorbar(nonlin.y - mean(nonlin.y), nonlin.x, nonlin.e, 'k');

        
        linkaxes([h1, h2], 'x');
        linkaxes([h2, h3], 'y');
        
        % shuffled control for sta
        n_repeats = 1000;
        cc_rand = zeros(n_repeats,1);
        sta_shuffle = zeros(n_repeats, 1);
        
        for i=1:n_repeats
            
            sta_shuffle = Utilities.Shuffle(x);
            
            cc_rand(i) = corr(sta_shuffle',sorted_resp);
            %mag=max(abs(gen2));
        end

        cc = corr(x',sorted_resp); %
        p = sum(cc_rand > cc)/length(cc_rand);
        
        % significance of projection
        h4 = subplot(3, 3, 1);
        h = histogram(cc_rand,0:0.01:1);
        h.FaceColor = [1 1 1 ];
        hold on;
        plot([cc cc],[0 50], 'LineWidth', 2, 'Color', 'r');
        text(.2, 40, num2str(p,'p = %.3f'))
        
        if greyImsOnly
            filename = [a_sessPath filesep  '2DGrid_' patID '_Cell_' num2str(aft_matches(cellIndex)) '_B&WIms'];  
        else
            filename = [a_sessPath filesep  '2DGrid_' patID '_Cell_' num2str(aft_matches(cellIndex)) '_ColoredIms'];   
        end
        
        print(f, filename, '-dpng', '-r0')

        close all
    end
    
end

%% helpers
function nonlin = compute_binned_average(rlin, resp, nbin, least_samples)


%% flexable:
% set an initial bin, combine near bins with too few samples
if nargin<3
    nbin = 8;
end
if nargin<4
    least_samples = 20;
end

mi=min(rlin);
ma=max(rlin);
edge=mi:(ma-mi)/nbin:ma; % min of proj to max of proj in steps of (ma-mi)/nbin

x=zeros(nbin,1);
y=zeros(nbin,1);
e=zeros(nbin,1);
ns=zeros(nbin,1); % number of samples

c = 0; % count

for i=1:nbin
    ind=rlin>edge(i) & rlin<=edge(i+1);
    ns(i)= sum(ind); %num of spikes between edges
end

edge_select = true(nbin+1,1);
for i=1:nbin
    if ns(i) < least_samples
        if i<nbin
            ns(i+1)=ns(i+1)+ns(i);
            edge_select(i+1) = false;
        elseif i==nbin
            ind = find(edge_select(1:end-1),1,'last');  % find last bin edge and combine the last bin to previous one
            edge_select(ind) = false;
        end
    end
end

edge = edge(edge_select); % deselect some dividing edges
nbin = length(edge)-1;
x=zeros(nbin,1);
y=zeros(nbin,1);
e=zeros(nbin,1);
ns=zeros(nbin,1);
for i=1:nbin
    ind=rlin>edge(i) & rlin<=edge(i+1); % images with sta_proj in this range
    x(i)=mean(rlin(ind));
    y(i)=mean(resp(ind));
    e(i)=std(resp(ind))/sqrt(sum(ind));
    ns(i)= sum(ind);
end

nonlin.x=x;
nonlin.y=y;
nonlin.e=e;
nonlin.ns=ns;
nonlin.rlin_std =std(rlin);
end
