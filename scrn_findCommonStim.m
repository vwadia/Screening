
%% Script to find common stimuli between 1593 object screening sessions
%% Set Paths

% task = 'Object_Screening';
task = 'Recall_Task';

pt_ID = 'P71CS';
% sess_ID = 'LargeObjectScreening_Session_1_20201123';
% sess_ID = 'FastObjectScreening_Session_1_20201125';
sess_ID = 'ObjectScreening_Session_1_20201121';

% pt_ID = 'P73CS';

[atCedars, basePath, ~] = Code.setPaths(task, pt_ID, sess_ID);


%% Read in images

pathStimuli = [basePath filesep 'stimuliUsed'];


imDir = dir(fullfile(pathStimuli));
imDir = imDir(~ismember({imDir.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));

for i = 1:length(imDir)
    name = imDir(i).name;
    imNames{i, 3} = str2num(name(1:end-4)); 
end

%%
% arr1 = imNames(:, 1);
% arr1 = arr1(1:434);
% 
% arr2 = imNames(:, 2);
% arr2 = arr2(1:500);
% 
% arr3 = imNames(:, 3);
% arr3 = arr3(1:200);

%%
overlap = intersect(cell2mat(imNames(:, 1)), cell2mat(imNames(:, 2)));