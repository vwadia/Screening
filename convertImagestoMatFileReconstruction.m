
% reading in images and then comverting them to the correct 4D mat file for
% the image reconstruction script

% it also serves to convert the generated images into a format (correct
% size and filename etc.) for synthetic screening

% 4D mat file's dimensions are #im x 3 x H x W - that's what pinglei used

% vwadia July2022

setDiskPaths

% path to images
% pathStimuli = 'G:\SUAnalysis\Object_Screening\P79CS\FingerprintScreening_Session_1_20220331\stimuliUsed';
% pathStimuli = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'FingerprintScreening_Session_1_20221026' filesep 'stimuliUsed'];
% pathStimuli = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030' filesep 'genStimTif_color_2'];
% pS2 = [diskPath filesep 'Object_Screening' filesep 'P81CS' filesep 'ClosedLoopScreening_Session_1_20221030' filesep 'genStimTif_color_2'];


pathStimuli = [diskPath filesep 'Object_Screening' filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115' filesep 'genStim'];
pS2 = [diskPath filesep 'Object_Screening' filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115' filesep 'effectiveStim_Synth'];

if ~exist(pS2)
    mkdir(pS2)
end

convertToTif = 1;
removeColor = 1;

num_steps = 10; % in each direction

if convertToTif

    imDir = Utilities.readInFiles(pathStimuli);
    
    [~, natIdx] = natsortfiles({imDir.name});
    imDir = imDir (natIdx); % don't really need this - but be aware
%     keyboard 
    
arrangement = [20:-1:16 11:1:15 10:-1:6 1:1:5]; % for pref and ortho
arr2 = [120:-1:116 110:115 109:-1:105 99:104 98:-1:94 88:93 87:-1:83 77:82 ...
    76:-1:72 66:71 65:-1:61 55:60 54:-1:50 44:49 43:-1:39 33:38 32:-1:28 22:27 ...
    21:-1:17 11:16 10:-1:6 1:5];
arr2_add = [1:55 111:120 100:110 89:99 78:88 67:77 56:66]; 
% arrange imdir
cell_one = imDir(1:20);
cell_one = cell_one(arrangement);

cell_two = imDir(21:40);
cell_two = cell_two(arrangement);

cell_three_1 = imDir(41:60);
cell_three_1 = cell_three_1(arrangement);
cell_three_2 = imDir(61:end);
cell_three_2 = cell_three_2(arr2);
cell_three_2 = cell_three_2(arr2_add);

imDir = [cell_one; cell_two; cell_three_1; cell_three_2];

    for i = 1:length(imDir)
        
% old way ---------------
%         name = imDir(i).name;
%         strpos = strfind(name, '.');
%         dashpos = strfind(name, '_');
%         better_name = name(9:strpos-1);
%         
%         if strcmp(name(dashpos+1), 'm')
%             
%             if (strpos - dashpos == 7)               
%                 stepnum = name(strpos-1);                
%             elseif (strpos - dashpos == 8)               
%                 stepnum = name(strpos-2:strpos-1);               
%             end
%             
%             realStepNum = ((num_steps + 1) - str2num(stepnum));
% 
%         elseif  strcmp(name(dashpos+1), 'p')
%             
%             if (strpos - dashpos == 6)
%                 stepnum = name(strpos-1);
%             elseif (strpos - dashpos == 7)
%                 stepnum = name(strpos-2:strpos-1);
%             end 
%             
%             realStepNum = num_steps + str2num(stepnum);
%             
%         end


        % new way --------------------
        name = imDir(i).name;
        dashpos = strfind(name, '_');
        dashpos = dashpos(end);
        n1 = name(1:dashpos-1);
        n2 = name(dashpos+1:end);
        
        
        
        
        I = imread([imDir(i).folder filesep imDir(i).name]);
        I = imresize(I, [224 224]);
        if removeColor
            
            fnum = sprintf('%03d', i);
            better_name = [n1 '_' fnum '_' n2];
            bestname = better_name;% ['cell' better_name(1:strfind(better_name, '_')-1) '_' fnum better_name(strfind(better_name, '_'):end)];
            
            
            imwrite(rgb2gray(I), [pS2 filesep [bestname '.tif']])
        else
            fnum = sprintf('%03d', length(imDir)+i);
            better_name = [n1 '_' fnum '_' n2];
            bestname = better_name;% ['cell' better_name(1:strfind(better_name, '_')-1) '_' fnum better_name(strfind(better_name, '_'):end)];
            imwrite(I, [pS2 filesep [bestname '.tif']])
        end
    end

end

%% name fixing 
setDiskPaths
pS2 = [diskPath filesep 'Object_Screening' filesep 'P82CS' filesep 'ClosedLoopScreening_Session_1_20230115' filesep 'effectiveStim_Synth'];

oldDir = pS2;
% oldDir = 'G:\SUAnalysis\Object_Screening\P82CS\ClosedLoopScreening_Session_1_20230115\All_images\sta_grid_904';

imDir = Utilities.readInFiles(oldDir);
[~, natIdx] = natsortfiles({imDir.name});
imDir = imDir (natIdx); % don't really need this - but be aware
newDir = [oldDir filesep 'NoSpaces'];
if ~exist(newDir)
    mkdir(newDir)
end

% getting rid of spaces in coordinate names - if Python is right won't need
% this at all
for i = 1:length(imDir)
    
    I = imread([imDir(i).folder filesep imDir(i).name]);
    curName = imDir(i).name;
    dashpos = strfind(curName, ',');
    n1 = curName(1:dashpos);
    n2 = curName(dashpos+1:end);
    if strcmp(n2(1), ' ')
        newName = [n1 n2(2:end)];
    else
        newName = [n1 n2];
        assert(isequal(curName, newName));
    end
    imwrite(I, [newDir filesep newName])
    
end


% adding 'grid' names to grid images
% for i = 1:length(imDir)
%     
%     I = imread([imDir(i).folder filesep imDir(i).name]);
%     curName = imDir(i).name;
%     dashpos = strfind(name, '_'); % NOTE THIS ERROR - ARE THE NAMES CORRECT? yes got lucky
%     n1 = curName(1:dashpos-1);
%     n2 = curName(dashpos+1:end);
%     
%     newName = [n1 '_Grid_' n2];
%     imwrite(I, [newDir filesep newName])
%     
% end




%%
% images need to have 3 channels
threeD = 1; 

% essentially just a vector from 1 -> # images
imageIDs = [1:500]';

grayImages = Utilities.getImageDescriptions(pathStimuli, 224, imageIDs, threeD);

imall = permute(grayImages, [4 3 1 2]);

save('varunImgs', 'imall');

%% display check - to make sure dimension permutation didn't mess with anything 

imgNum = 55;
sampleIm = squeeze(imall(imgNum, :, :, :));

s_Im = squeeze(sampleIm(1, :, :));

figure; imshow(s_Im, [])

%% double display check -  with pinglei's images. permute them back

load('G:\SUAnalysis\ObjectSpace\ImageReconstructionCode_Python\imalltest2.mat')

imall_perm = permute(imall, [3 4 2 1]);

imgNum = 25;
sampleIm = squeeze(imall(imgNum, :, :, :));

s_Im = squeeze(sampleIm(1, :, :));

figure; imshow(s_Im, [])