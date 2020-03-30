%% Analyze real experiment data

%% This should be run on atlas, no local copy
clear all;close all;clc

%% Here it is  
cd /home/range1-raid2/sgwarren-data/Deviant/Blofeld/Data/Deviant/B % monkey B
% cd /home/range1-raid2/sgwarren-data/Deviant/Kananga/Data/Deviant/K monkey

on = load('DeviantOn_AlignTo_stimulusOnTime.mat');

cd ~/Dropbox/stonesync/19attentionprobV1optimaging/opticalattention/
onImg = ImageHelper.convertSparseToFull(on.S, on.IX, on.V);

nOnTrials = size(onImg, 4);

%% set parameters
timeWindow = 5:20; % B, ON:5:20; OFF,3:15;
stimOnset = 11; % B, ON:11,OFF,6;
imgMask = []; % image Mask, [xrange; yrange];
wantdemean = 1; % whether to demean the data 

%% analyze behavioral data
correctIdx = zeros(1,nOnTrials);
deviantPos = zeros(1,nOnTrials);
for i = 1:nOnTrials
    correctIdx(i) = on.T(i).behavior; % 0,correct; 1, false positive,2, false negative
    deviantPos(i) = on.T(i).trialDescription.deviantPosition;
end
% group the trial
[deviantTrialCnt, correctTrialCnt, hitRate] = deal(zeros(1,9));
for i=1:9
    deviantTrialCnt(i) = sum(deviantPos==(i-1));
    correctTrialCnt(i) = sum(deviantPos==(i-1) & correctIdx==0);
    hitRate(i) = correctTrialCnt(i) / deviantTrialCnt(i);
end

% Visualize the position sensitivity
close all;
h = figure;
set(h, 'Position', [0 0 700 600]);
% plot dots
hitRate = hitRate([9 3 6 8 2 5 7 1 4]); % rearange the position for visualization purpose
dotX = [1 1 1 2 2 2 3 3 3];
dotY = [3 2 1 3 2 1 3 2 1];
dotText = {'9', '3', '6', '8', '2', '5', '7', '1', '4'};

scatter(dotX, dotY, 6000, hitRate,'filled', 'MarkerEdgeColor','k');
text(dotX,dotY,dotText);
xlim([0 4]);ylim([0 4]);
cmin = 0.4;
cmax = 0.7;
colormap(h,jet(256)); caxis([cmin, cmax]); colorbar();
axis off;
set(gca,'View',[45 90]);
title('Positional senstivity attend in');

%savefig('B_posisensi_ON.fig');


%% Read out the trials
cueIdx = zeros(1,nOnTrials);
for i = 1:nOnTrials
    cueIdx(i) = on.T(i).trialDescription.explicitCue;
end

% Now group image into cells based on stimulusIndex
imgON = cell(1, length(unique(cueIdx)));
for i = [1 2]
    imgON{i} = onImg(:,:,:, cueIdx==(i-1));
end
%imgON{1}, attendin; imgON{2}, attendout;
%
imgON = cellfun(@(x) x(:,:,timeWindow,:), imgON, 'UniformOutput', 0); % extract time range
imgON = cellfun(@(x) nanmean(x, 4), imgON, 'UniformOutput', 0); % average across trials
% reshape the image to length*width x time_points
imgON = cellfun(@(x) reshape(x,[316*316 length(timeWindow)]), imgON, 'UniformOutput',0); 
validMask = cellfun(@(x) ~any(isnan(x),2), imgON, 'UniformOutput',0); % create a valid value mask, necessary to use it later
imgON = cellfun(@(x,y) x(y,:)', imgON, validMask, 'UniformOutput',0); % time X pixel*pixel in each cell

if wantdemean
    imgON = cellfun(@(x) x-repmat(mean(x,2),1,size(x,2)), imgON, 'UniformOutput',0);
end

%% load the PCs obtained from localizer experiment, and the PC obtained from PC on main experiment
localizerData = load('localizer_PCA_B.mat');
localizerMaps = localizerData.valueImg_mat_demean;
localizerMask = localizerData.nanMask;

PC1 = load('mainExp_PCA_B_ON.mat','PC1');
PC1 = PC1.PC1;
% Note that the individual maps, are conditions [9 3 6 8 2 5 7 1 4] for
% visualization purpose

%% project main data on the PC obtained from localizer data
valueLoad = cellfun(@(x) PC1'*x, imgON, 'UniformOutput',0); 
% valueLoad is a 1 x 9 cell, each element is a pixels X 1 vector

valueImg = repmat({NaN(1,316*316)},[1,2]);
valueImg = cellfun(@(x) reshape(x, [316 316]), valueImg, 'UniformOutput',0);
for i=1:2; valueImg{i}(validMask{i})=valueLoad{i}; end
% concatenate all images
valueImg_mat = cat(3, valueImg{:}); % concatenate images
%valueImg_mat_demean = valueImg_mat-repmat(nanmean(valueImg_mat,3),[1 1 2]);

% visualized denoised data
h = figure;
set(h,'Position',[0 0 800 300]);
imagesc(makeimagestack(valueImg_mat)); colorbar(); caxis([-2, 2]);
% we need to rerange the order of the figure
savefig(h, 'B_mainExp_denoisedmaindata.fig');


%% perform regression analysis

% Combine localizer masks and main exp masks
tmp = cat(2,validMask{:}, localizerMask{:});
validMaskAll = all(tmp,2); % validMaskAll is the mask for valid voxels for both localizer and main exp data.

% Extract individual maps data
localizerMaps = reshape(localizerMaps, [316 * 316, 9]);
localizerMapValid = localizerMaps(validMaskAll,:);
localizerMapValid(:,end+1) = 1;  % add a constant regressor here

% extract main data
valueImg_mat = reshape(valueImg_mat,[316 * 316, 2]);
valueImg_mat = valueImg_mat(validMaskAll, :);

% attend in
betaAttIn = valueImg_mat(:,1)'/localizerMapValid';
betaAttOut = [0.0836    0.2231   -0.0757    0.0439   -0.0352         0    0.3075    0.0273   -0.0515    0.0629];

% attend in
betaAttOut = valueImg_mat(:,2)'/localizerMapValid';
betaAttIn = [0.1948    0.3195   -0.0241    0.0953    0.0303         0    0.3069    0.1025   -0.0396    0.0387];
% 
betaAttInOutDiff = betaAttIn - betaAttOut;
% condition is [9 3 6 8 2 5 7 1 4]


% Visualize the difference map
dotX = [1 1 1 2 2 2 3 3 3];
dotY = [3 2 1 3 2 1 3 2 1];
dotText = {'9', '3', '6', '8', '2', '5', '7', '1', '4'};

close all;
h = figure;
set(h, 'Position', [0 0 700 600]);
% plot background
cmin = -0.12;
cmax = 0.12;
cmap=colormap(h,jet(256)); caxis([cmin, cmax]);
valueBg = betaAttInOutDiff(end);
index = fix((valueBg-cmin)/(cmax-cmin)*size(cmap,1))+1;
cBg = cmap(index,:);
patch([0 4 4 0], [0 0 4 4],cBg); 
colorbar(); hold on;

% plot dots 
scatter(dotX, dotY, 6000, betaAttInOutDiff(1:end-1)','filled', 'MarkerEdgeColor','k');
text(dotX,dotY,dotText);
xlim([0 4]);ylim([0 4]);
axis off;
set(gca,'View',[45 90]);
title('attend in - out');

savefig('B_attend_betadiff_ON.fig');


