
setDiskPaths 

% OG images
load([diskPath filesep 'ObjectSpace' filesep 'ImageReconstructionCode_Python' filesep 'varunImgs.mat']);

% generated images
load([diskPath filesep 'ObjectSpace' filesep 'ImageReconstructionCode_Python' filesep 'varunOutput1.mat']);



%% noise injection
% n_pixels = numel(imall(1, 1, :, :)); % should give you 227*227 = 51529
% 
% percent = 10;
% 
% n_zeros = floor((20/100)*n_pixels); % how many 0's you want this is 10% for now
% 
% noisyIms = zeros(size(imall));
% 
% tic % this will time the section of code within 'tic' and 'toc'
% for im = 1:size(imall, 1) % number of images    
%     
%     sample_im = squeeze(imall(im, :, :, :));
%     
%     sample_im_wNoise = sample_im;
%     
%      % randomly choose a fixed number of indices (n_zeros) to 0 
%     sample_im_wNoise(randperm(n_pixels, n_zeros)) = 0;
%     
%     % store each altered image
%     noisyIms(im, :, :, :) = sample_im_wNoise;
% end
% toc

% save
% save(['varunImgs_' num2str(percent) 'percent_noise'], 'noisyIms')

%% 

sim = [];
sim2 = [];

for im = 1:size(imall, 1)
    
    A = squeeze(noisyIms(im, :, :, :));
    REF = squeeze(imall(im, :, :, :));
    sim(im) = ssim(REF, A);
    sim2(im) = ssim(A, REF);
end
figure; hold on; scatter(1:500, sim); scatter(1:500, sim2)
title('noise 20%') 