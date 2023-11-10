% script to analyze catch question responses for all sessions *fuck me*
% First:
%     Create lookup table (10 possible questions - I think)
%         As a test - first compile questions, then session list
%         loop through sessions, load taskstruct and match questions asked with the q's
%         If any have questions that don't match the 10 keep track of those sessions/questions
%
%     Add those questions to the lookup table
%     (All ims in same category should have roughly the same answers to quesitons)
%
% Second:
%     Write code to loop through sessions, load taskStruct, compare all images and find idx #s for them
%     use lookup table to get correct/wrong answers and then percentage of
%     total.
%     load in events too and record response times? (Resp - im off, I don't send TTL for catch q. Change this in future iterations)

setDiskPaths

% task = 'Object_Screening';
% task = 'Recall_Task';
task = 'ReScreen_Recall';

% current set of catch questions
options{1} = 'Was the last image of an animal?';
options{2} = 'Did the last image contain a plant/fruit/vegetable?';
options{3} = 'Could object in the last image fit in your palm?';
options{4} = 'Was there a person in the last image?';
options{5} = 'Was the last image of a building or vehicle?';
options{6} = 'Was the object in the last image a letter or a number?';
options{7} = 'Was the object in the last image a living thing?';
options{8} = 'Was the last image a cartoon (as opposed to a photograph)?';
options{9} = 'Was the last image black and white (as opposed to color)?';
options{10} = 'Did the last image contain a musical instrument?';
options{11} = 'Could the object in the last image fly?';

options{12} = '¿Fue la última imagen un animal?';
options{13} = '¿La última imagen contenía una planta/fruta/verdura?';
options{14} = '¿Podría el objeto de la última imagen caber la palma de tu mano?';
options{15} = '¿Había una persona en la última imagen?';
options{16} = '¿Fue la última imagen de un edificio o un vehículo?';
options{17} = '¿El objeto de la última imagen era una letra o un número?';
options{18} = '¿Era el objeto de la última imagen un ser vivo?';
options{19} = '¿Fue la última imagen una caricatura (a diferencia de una fotografía)?';
options{20} = '¿La última imagen era en blanco y negro (en lugar de color)?';
options{21} = '¿La última imagen contenía un instrumento musical?';
options{22} = '¿Podría volar el objeto de la última imagen?';

options = options';


imageIDs = [1:500]';


%% build session list


[sessID, tStrct, ~] = Utilities.sessionListAllTasks(task, true, false);

%% Step 1 check for extra catch questions and add them

new_qs = {};

for ss = 1:length(sessID)
    
    load([diskPath filesep sessID{ss} filesep tStrct{ss}])
    
    % catchResponses
    % images_questioned
    % questions_asked
    new_qs = cat(1, new_qs, setdiff(questions_asked', options));
    
    
end
new_qs = unique(new_qs);  % there are repetitions

% if strcmp(task, 'Object_Screening')
%     catchQs = cat(1, options', new_qs);
% else
catchQs = cat(1, options, new_qs);
% end

%% create lookup table (partially manual)

% lookUpTable_values = false(length(catchQs), length(imageIDs));
%
% % read in Images - from any session
% ims = Utilities.readInFiles([diskPath filesep 'Object_Screening' filesep '500Stimuli']);
% images = {};
% for i = 1:length(ims)
%
%     % match the way they are in taskStruct
%     images{i} = imread([ims(i).folder filesep ims(i).name]);
%
% end
%
% % set up category orders
% faceInds = 134:210;
% objInds = [85:133 236:255 283:289 291:356 409:500]; % chnged to include 290 in text vwadia march 2022
% textInds = [264:282 290 400:408];
% vegInds = [211:235 357:399];
% animInds = [1:84 256:263];
%
% catOrder = zeros(length(order), 1);
% catOrder(ismember(order, faceInds)) = 1;
% catOrder(ismember(order, textInds)) = 2;
% catOrder(ismember(order, vegInds)) = 3;
% catOrder(ismember(order, animInds)) = 4;
% catOrder(ismember(order, objInds)) = 5;
%
% %% actually fill out look up table
%
% % category - english
% lookUpTable_values(1, animInds) = true; % animals
% lookUpTable_values(2, vegInds) = true; % plant or fruit
% lookUpTable_values(4, faceInds) = true; % faces
% lookUpTable_values(6, textInds) = true; % letter or number
% lookUpTable_values(7, [faceInds vegInds animInds]) = true; % living things
% lookUpTable_values(9, :) = true; % black or white photo
%
% % category - Spanish
% lookUpTable_values(12, animInds) = true; % black or white photo
% lookUpTable_values(13, vegInds) = true; % plant or fruit
% lookUpTable_values(15, faceInds) = true; % faces
% lookUpTable_values(17, textInds) = true; % letter or number
% lookUpTable_values(18, [faceInds vegInds animInds]) = true; % living things
% lookUpTable_values(20, :) = true; % black or white photo
%
% lookUpTable_values(24, :) = false; % scrambled background
%
% lookUpTable_values(33, animInds) = true; % animals
% lookUpTable_values(35, faceInds) = true; % faces
% lookUpTable_values(36, :) = false; % throw at me
%
% % manually check
%
% % fit in palm
% lookUpTable_values(3, [11 25 28 29 34 45 46 59 60 68 71:73 87:89 103 107:110 ...
%     112 113 115:117 119 122 124 125 126 129 132 133 211:225 227:235 255:263 286 ...
%     292 297 298 301 303 304 310 314 317:319 322 324:329 332 335 336 339 340 ...
%     344 357 358 361 362 369 375 376 382 383 384 385 389 391 393 394 396 399 ...
%     412 414 417 420:426 429 430 433 435 442 451 452 456 460 469 470 472:475 477 478 ...
%     479]) = true;
% lookUpTable_values(14, [11 25 28 29 34 45 46 59 60 68 71:73 87:89 103 107:110 ...
%     112 113 115:117 119 122 124 125 126 129 132 133 211:225 227:235 255:263 286 ...
%     292 297 298 301 303 304 310 314 317:319 322 324:329 332 335 336 339 340 ...
%     344 357 358 361 362 369 375 376 382 383 384 385 389 391 393 394 396 399 ...
%     412 414 417 420:426 429 430 433 435 442 451 452 456 460 469 470 472:475 477 478 ...
%     479]) = true;
%
% % building or Vehicle
% lookUpTable_values(5, [249 252 253 461 462 480:485 488 489 494:496 500]) = true;
% lookUpTable_values(16, [249 252 253 461 462 480:485 488 489 494:496 500]) = true;
%
% % cartoon or photo
% lookUpTable_values(8, [3 4 13 20 21 32 55 84 108 249 252 textInds 287 289 290 ...
%     295 296 302 305 309 315 316 337 338 356 359 360 364:368 373 375:376 378 ...
%     380 381 397 399 412 416 417 433 439 440 480 496]) = true;
% lookUpTable_values(19, [3 4 13 20 21 32 55 84 108 249 252 textInds 287 289 290 ...
%     295 296 302 305 309 315 316 337 338 356 359 360 364:368 373 375:376 378 ...
%     380 381 397 399 412 416 417 433 439 440 480 496]) = true;
% lookUpTable_values(25, [3 4 13 20 21 32 55 84 108 249 252 textInds 287 289 290 ...
%     295 296 302 305 309 315 316 337 338 356 359 360 364:368 373 375:376 378 ...
%     380 381 397 399 412 416 417 433 439 440 480 496]) = true;
% lookUpTable_values(32, [3 4 13 20 21 32 55 84 108 249 252 textInds 287 289 290 ...
%     295 296 302 305 309 315 316 337 338 356 359 360 364:368 373 375:376 378 ...
%     380 381 397 399 412 416 417 433 439 440 480 496]) = true;
% lookUpTable_values(37, [3 4 13 20 21 32 55 84 108 249 252 textInds 287 289 290 ...
%     295 296 302 305 309 315 316 337 338 356 359 360 364:368 373 375:376 378 ...
%     380 381 397 399 412 416 417 433 439 440 480 496]) = true;
%
% % musical instrument
% lookUpTable_values(10, [ 122 341:352 354 443 445]) = true;
% lookUpTable_values(21, [ 122 341:352 354 443 445]) = true;
% lookUpTable_values(28, [ 122 341:352 354 443 445]) = true;
%
% % fly
% lookUpTable_values(11, [4 6 15 24 26 32 54:56 59 60 71 72 74 108 256:259 ...
%     485 494:496]) = true;
% lookUpTable_values(22, [4 6 15 24 26 32 54:56 59 60 71 72 74 108 256:259 ...
%     485 494:496]) = true;
%
% % plant only
% lookUpTable_values(23, [357 358 370:372 379 384 386:388]) = true;
% lookUpTable_values(38, [357 358 370:372 379 384 386:388]) = true;
%
% % fruit only
% lookUpTable_values(29, [211:235 361 362]) = true;
% lookUpTable_values(34, [211:235 361 362]) = true;
%
%
% %% shape judgements using PC1PC2 spread because I am lazy/it is better ground truth
% % quadrants 2 and 3 = stubby
% % quadrants 1 and 4 = spiky
% % x value is what matters
%
% % %       |
% % %    2  |  1
% % % ------+------
% % %    3  |  4
% % %       |
%
%
% load([diskPath filesep 'ObjectSpace' filesep '500Stimuli' filesep 'params_Alexnet_fc6_500Stimuli.mat']);
% if exist('params')
%     score = params;
% end
% PC12 = score(:, 1:2);
% H1 = []; H2 = [];
%
% for ii = 1:length(PC12(:, 1))
%
%     if PC12(ii, 1) > 0 % spiky
%         H1 = vertcat(H1, ii);
%     elseif PC12(ii, 1) < 0 % stubby
%         H2 = vertcat(H2, ii);
%     end
% end
%
%
% % labels
% PC12(H1, 3) = 1;
% PC12(H2, 3) = 2;
%
%
% % pointy or not...fuck me use PC1:PC2 rep to answer these
% lookUpTable_values(30, [find(PC12(:, 3) == 1)]) = true;
%
% % round or not...what a dumb question
% lookUpTable_values(31, [find(PC12(:, 3) == 2)]) = true;
%
% % save([diskPath filesep 'Object_Screening' filesep 'CatchQLookUpTable'], 'lookUpTable_values')

%% Now the fun. Load in the images compare thenm and actually find the naswwers

setDiskPaths
load([diskPath filesep 'Object_Screening' filesep 'CatchQLookUpTable']);

sessCatchPerf = nan(length(sessID), 1); % % of correct catch q answers for each session
% sessCatchPerf = []; % % of correct catch q answers for each session

% read in Images - from any session
ims = Utilities.readInFiles([diskPath filesep 'Object_Screening' filesep '500Stimuli']);
images = {};
for i = 1:length(ims)
    
    % match the way they are in taskStruct
    images{i} = imread([ims(i).folder filesep ims(i).name]);
    
end

% for each session
for ss = 1:length(sessID)
    
    % load taskStrcut
    load([diskPath filesep sessID{ss} filesep tStrct{ss}])
    
    
    % load in images for that session- have to do it this way
    % for close loop sessions with synth ims
    if strcmp(task, 'ReScreen_Recall')
        if ss == 9 || ss == 10
            ims = Utilities.readInFiles([diskPath filesep sessID{ss} filesep 'stimuliUsed']);
            images = {};
            
            for i = 1:length(ims)
                % match the way they are in taskStruct
                images{i} = imread([ims(i).folder filesep ims(i).name]);
                
            end
        end
    end
    
    if strcmp(task, 'Object_Screening') && ss == 1 % the images didn't get saved here
        continue
    else
        % some early sessions got fucked
        catchAns = [];
        givenAns = []; corrAns = [];
        ctr = 1;
        if iscell(images_questioned{1})
            
            for im_q = 1:length(images_questioned)
                
                idx_im = cellfun(@(x) isequal(images_questioned{im_q}{1}, x), images, 'UniformOutput', false);
                idx_im = find(cell2mat(idx_im));
                
                if idx_im <=500
                    
                    % question asked
                    idx_q = cellfun(@(x) isequal(questions_asked{im_q}, x), catchQs, 'UniformOutput', false);
                    idx_q = find(cell2mat(idx_q));
                    
                    % correct answer
                    corrAns(ctr) = lookUpTable_values(idx_q, idx_im);
                    
                    
                    if strcmp(catchResponses{im_q}, 'False')
                        givenAns(ctr) = 0;
                    elseif strcmp(catchResponses{im_q}, 'True')
                        givenAns(ctr) = 1;
                    end
                    
                    catchAns(ctr) = isequal(corrAns(ctr), givenAns(ctr));%lookUpTable_values(idx_q, idx_im);
                    ctr = ctr + 1;
                end
            end
            
            
        else %do it by index number (for P73)
            im = cell2mat(images_questioned)';
            im = mod(im, 1593);
            imNames = cat(1, ims(:).name); imNames = str2num(imNames(:, 1:end-4));
            
            for im_q = 1:length(im)
                idx_im = im(im_q);
                if idx_im <=500
                    
                    % question asked
                    idx_q = cellfun(@(x) isequal(questions_asked{im_q}, x), catchQs, 'UniformOutput', false);
                    idx_q = find(cell2mat(idx_q));
                    
                    % correct answer
                    corrAns(ctr) = lookUpTable_values(idx_q, idx_im);
                    
                    
                    if strcmp(catchResponses{im_q}, 'False')
                        givenAns(ctr) = 0;
                    elseif strcmp(catchResponses{im_q}, 'True')
                        givenAns(ctr) = 1;
                    end
                    
                    catchAns(ctr) = isequal(corrAns(ctr), givenAns(ctr));%lookUpTable_values(idx_q, idx_im);
                    ctr = ctr + 1;
                end
            end
        end
        % get sess average
        sessCatchPerf(ss) = sum(catchAns)/length(catchAns);
    end
end

% I inverted the colors on the button box for these sessions
% P80 fingerprint session 1
% P87 fingerprint screen
if strcmp(task, 'Object_Screening')
    sessCatchPerf(15) = 1 - sessCatchPerf(15);
    sessCatchPerf(30) = 1 - sessCatchPerf(30);
end


% have to run for object screening and rescreen recall and concatenate both
% sessCatchPerf_All = [];
% sessCatchPerf_All = cat(1, sessCatchPerf_All, sessCatchPerf);
% clearvars -except sessCatchPerf_All

% save([diskPath filesep 'Object_Screening' filesep 'PerSessionCatchQResponses'], 'sessCatchPerf_All');

%% plotting

setDiskPaths

% for grouping - eg. sessions 4 - 8 were the same patient
patIdx = [1 2 3 4 4 4 4 5 6 6 7 7 7 7 8 8 8 8 9 9 9 10 10 11 11 11 12 12 13 14 15 4 4 4 7 7 7 8 8 9 10 11 11 12 12];
allPatsTogether = false;

load([diskPath filesep 'Object_Screening' filesep 'PerSessionCatchQResponses'])
load([diskPath filesep 'Object_Screening' filesep 'CatchQLookUpTable']);

err = nanstd(sessCatchPerf_All)./sqrt(sum(~isnan(sessCatchPerf_All)));
X = nanmedian(sessCatchPerf_All);

f = figure;
set(f, 'Position', [680 550 600 450]) % so the xlabels aren't tilted 
hold on
title('Catch question responses')
filename = [diskPath filesep 'Object_Screening' filesep 'forPaper' filesep 'CatchQ_Summary'];

if allPatsTogether
    
    % bar graph
    b = bar(1, X, 0.06);
    errorbar(1, X, err, 'k', 'LineStyle', 'none')
    b.FaceColor = [0.8500 0.45 0.05];
    b.FaceAlpha = 0.8;
    b.EdgeAlpha = 0.1;
    
    % scatter
    scatter(ones(1, length(sessCatchPerf_All)), sessCatchPerf_All, 20, [0 0.2 0], 'filled');
    
    % boxplot
    boxplot(sessCatchPerf_All, 'BoxStyle','outline', 'Widths', 0.015)
    
    % text
    str = {['Median = ' num2str(X*1e2,' %.2f') '%']};
    text(1.02, 0.89, str, 'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold')
    
    str2 = {['Sessions = ' num2str(length(sessCatchPerf_All))]};
    text(1.02, 0.92, str2, 'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold')
    
    xlim([0.9 1.1])
    ylim([0.5 1])
    
else
    % remove the 2 pts with sessions where catch images weren't saved
    badIdx = find(isnan(sessCatchPerf_All));
    patIdx = patIdx(badIdx+1:end)-length(badIdx);
    sessCatchPerf_All = sessCatchPerf_All(badIdx+1:end);
    % scatter
    scatter(patIdx, sessCatchPerf_All, 20, [0 0.2 0], 'filled');
    
    % boxplot
    boxplot(sessCatchPerf_All, patIdx,  'BoxStyle','outline');%, 'Widths', 0.015)
    
    % line
    yline(0.5, '--k', 'Linewidth', 2)
    
    % text
    str = {['Median = ' num2str(X*1e2,' %.2f') '%']};
    text(5.5, 0.35, str, 'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold')
    
    str2 = {['Sessions = ' num2str(length(sessCatchPerf_All)+2)]};
    text(5.5, 0.3, str2, 'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold')
    
    text(0.55, 0.53, 'Chance', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold')
    
    ylim([0 1])
    
    filename = [filename '_perPat'];
    xlabel('Patient ID')
    
end
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

ylabel('Proportion of correct answers')

print(f, filename, '-dpng', '-r300')
close all


