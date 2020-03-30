% PCA analysis 

clear all;close all;clc;

%% load data
%cd /home/range1-raid2/sgwarren-data/Deviant/Kananga/Data/Deviant/K  % go to the data folder
%datadir = '/home/range1-raid3/sgwarren/matlab/Deviant/Data/K';

cd /home/range1-raid2/sgwarren-data/Deviant/Blofeld/Data/Deviant/B
datadir = '/home/range1-raid2/sgwarren-data/Deviant/Blofeld/Data/Deviant/B';

% load the data
posidata = load('PositionTuning_AlignTo_stimulusOnTime.mat'); 


cd ~/Dropbox/stonesync/19attentionprobV1optimaging/opticalattention/
% convert images from sparse to full using scott's function
posiTuningImg = ImageHelper.convertSparseToFull(posidata.S,posidata.IX, posidata.V);
%posiTuningImg a 316 X 316 x 33(time points) x n(trials)

%% other setup
timeWindow = 5:20; 
% between all 1 x 33 time points, we only want to extract part of them
% Monkey B: 5-20; Monkey K: 3-15

stimOnset = 11; % stimulus onset frame
% Frame 11 stimulus onset for Monkey B; frame 6 stimulus onset for Monkey K;

% You may want to mask out some pixels, like "" 
%imgMask = [1 316;30 120]; % image Mask, [xrange; yrange];
imgMask = []; % image Mask, [xrange; yrange];

wantdemean = 1; % whether to demean the data for each condition before perform PCA

K = 2; % how many principle component you want to obtan

%% Preprocessing 1
nTrials = size(posiTuningImg, 4); % total number of trials

% extract condition index
allStimIndex = zeros(1,nTrials);
for i=1:nTrials
    allStimIndex(i) = posidata.T(i).trialDescription.deviantPosition;
end

% Group images into different conditions
img = cell(1,18); % we total have 18 conditions, 0~17
for i=1:18 % 18 positions
    img{i} = posiTuningImg(:, :, :, allStimIndex==(i-1));
end
% now img is a 1 x 18 cell, each element is a length x width x time_points
% x trial matrix

otherSideImg = img(10:18); % Condition 10-18 is from the other side

%===== note here======
%posiTuningImg = img(1:end); % We only need the first 9 positions
posiTuningImg = img(1:9); % We only need the first 9 positions
% Here we should decide whether we use all 18 conditions and only first 9
% conditions

%posiTuningImg = otherSideImg;
posiTuningImg = cellfun(@(x) permute(x,[1 2 4 3]), posiTuningImg, 'UniformOutput',0); % height X width X trial X time;
posiTuningImgCopy = posiTuningImg;

%% Preprocessing 2

% average across trials, here we use nanmean because there is a lot of nan
% pixels
posiTuningImg = cellfun(@(x) squeeze(nanmean(x, 3)), posiTuningImgCopy, 'UniformOutput',0); % average trials

% extract time window
posiTuningImg = cellfun(@(x) x(:, :, timeWindow), posiTuningImg, 'UniformOutput',0); % extract time window;

% mask out some pixels, set them to nan
if ~isempty(imgMask)
    for i=1:9, posiTuningImg{i}(imgMask(1,1):imgMask(1,2), imgMask(2,1):imgMask(2,2), :)=nan; end %
end

% reshape the image to length*width x time_points
posiTuningImg = cellfun(@(x) reshape(x,[316*316 length(timeWindow)]), posiTuningImg, 'UniformOutput',0); 

% Exclude NaN pixels and transpose
nanMask = cellfun(@(x) ~any(isnan(x),2), posiTuningImg, 'UniformOutput',0); % create a nan value mask, necessary to use it later
posiTuningImg = cellfun(@(x,y) x(y,:)', posiTuningImg, nanMask, 'UniformOutput',0); % time X pixel*pixel in each cell
%now each element is a time_points X pixel*pixel matrix

% In each condition, remove the mean activity, help remove effects of bubbles.
if wantdemean
    posiTuningImg = cellfun(@(x) x-repmat(mean(x,2),1,size(x,2)), posiTuningImg, 'UniformOutput',0);
end

%% --------------- do grand PCA --------------------- 
% We perform a PCA on the combined images from the average-trial images
% across all conditions

% combine all avg images
grandImg = cat(2,posiTuningImg{:});

% run singular vector decomposition
% the U is the eigen vectors sorted by its eigen values
[U,S,VT] = svd(grandImg,'econ');

% Let's visualize the first K components
close all;
h=figure;
set(h,'Position',[0 0 400 300]);
lh=myplot((timeWindow-stimOnset) * 0.2, U(:,1:K)');
straightline(0,'v','k'); % add the stimulus onset line
xlabel('SOA (sec)');ylabel('Signal change (%)');
set(lh(1),'Color','r');
set(lh(2),'Color','b');
legend(lh,{'PC1','PC2'}, 'Box','off');
set(gca,'Color','none');
%savefig(h,'B_localizer_PC_18cond.fig');
savefig(h,'B_localizer_PC_9cond.fig');

%% calculate the projected the average images for the 1st pc
compoGrand = -U(:,1); % We only project the first component

% posiTuningImg is a 1 x 9 cell, each element is a time_points x pixels
% matrix.
valueLoad = cellfun(@(x) compoGrand'*x, posiTuningImg, 'UniformOutput',0);
% valueLoad is a 1 x 9 cell, each element is a 1 X pixels vector

% But we have to add NaN pixels back to images such that we can fully construct the images 
valueImg = repmat({NaN(1,316*316)},[1,9]);
valueImg = cellfun(@(x) reshape(x, [316 316]), valueImg, 'UniformOutput',0);
for i=1:9; valueImg{i}(nanMask{i})=valueLoad{i}; end
% valueImg is a 1 x 9 cell, each element is a 316 x 316 matrix

% concatenate all iamges
valueImg_mat = cat(3, valueImg{[9 3 6 8 2 5 7 1 4]}); % concatenate images

% substract the mean image in each image
valueImg_mat_demean = valueImg_mat-repmat(nanmean(valueImg_mat,3),[1 1 9]);

% we crop the map and only keep the activation area
valueImg_mat_demean2 = valueImg_mat_demean(30:316-30,30:316-30,:);
%% Visualize the projection on the 1st component
close all;
h=figure;
set(h,'Position',[0 0 800 600]);
imagesc(makeimagestack(valueImg_mat_demean2)); % we need to rerange the order of the figure
set(gca,'visible','off');
caxis([-0.5, 0.5]); c = colorbar;
c.label.string='Signal change (%)';
savefig(h,'B_localizer_map_crop.fig');
%print(h,'-dpdf','-painters','-r300','vxssimu_estimation_calclfi1.pdf'); %save the figure to pdf


% ======= below is optinal ===========
% Demean the other side as a comparison
otherSideImg2 = cellfun(@(x) nanmean(nanmean(x,4),3), otherSideImg, 'UniformOutput', 0);
otherSideImg2 = nanmean(cat(3, otherSideImg2{:}),3);
valueImg_mat_demeanSide = valueImg_mat-otherSideImg2;
close all;
h = figure;
imagesc(makeimagestack(valueImg_mat_demean)); 
caxis([-0.5, 0.5]);colorbar();

%print(h,'-dpdf','-painters','-r300','vxssimu_estimation_calclfi1.pdf'); %save the figure to pdf


% %% project the single trial images for the 1st pc
% timeWindow2 = 2:11; % K, 2:11; B:4:12;
% compoGrand2 = compoGrand(timeWindow2); % further trancate time points
% 
% posiTuningImg = cellfun(@(x) x(:, :, :, timeWindow), posiTuningImgCopy, 'UniformOutput',0); % extract time window;
% if ~isempty(imgMask)
%     for i=1:9, posiTuningImg{i}(imgMask(1,1):imgMask(1,2), imgMask(2,1):imgMask(2,2), :,:)=nan; end % mask out the pixels
% end
% 
% % Reshape the matrix
% tmp = cellfun(@(x) reshape(x, [316*316*size(x,3), length(timeWindow)]), posiTuningImg, 'UniformOutput', 0);
% % Project in details
% posiTuningImg_proj_trial = cellfun(@(x) x(:,timeWindow2)*compoGrand2, tmp, 'UniformOutput', 0);
% % Reshape back
% posiTuningImg_proj_trial = cellfun(@(x) reshape(x, [316, 316, length(x)/(316*316)]), posiTuningImg_proj_trial, 'UniformOutput', 0);
% % Take the mean
% posiTuningImg_proj_trial_avg = cellfun(@(x) nanmean(x,3), posiTuningImg_proj_trial, 'UniformOutput', 0);
% % cat
% posiTuningImg_proj_trial_avg = cat(3, posiTuningImg_proj_trial_avg{:});
% imagesc(makeimagestack(posiTuningImg_proj_trial_avg));caxis([-1,1]);colorbar;

%%

