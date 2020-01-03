%% Analyze real experiment data

%% This should be run on atlas, no local copy
clear all;close all;clc

%% Here it is
cd /home/range1-raid2/sgwarren-data/Deviant/Blofeld/Data/Deviant/B % monkey B
% cd /home/range1-raid2/sgwarren-data/Deviant/Kananga/Data/Deviant/K monkey

on = load('DeviantOn_AlignTo_stimulusOnTime.mat');
onImg = ImageHelper.convertSparseToFull(on.S, on.IX, on.V);
nOnTrials = size(onImg, 4);


%% Read out the trials
deviantPos_ON = zeros(1,nOnTrials);
for i = 1:nOnTrials
    deviantPos_ON(i) = on.T(i).trialDescription.stimulusIndex;
end
% Now group image into cells based on stimulusIndex
stimulusCondON = unique(deviantPos_ON);
imgON = cell(1, length(stimulusCondON));
for i = 1:length(stimulusCondON)
    imgON{i} = onImg(:,:,:,deviantPos_ON==(i-2));
end
% We save the average image for each condition
allImg = cellfun(@(x) nanmean(nanmean(x,4),3), imgON, 'UniformOutput', 0);
allImg = cat(3, allImg{:});

close all;
h = figure;
set(h, 'Position', [0 0 800 600]);
imagesc(makeimagestack(allImg));colorbar(); caxis([-0.5, 0.5]);

%savefig(h, 'B_main_ON_avgmap.fig');


%% set PCA parameters
timeWindow = 5:20; % B, ON:5:20; OFF,3:15;
stimOnset = 11; % B, ON:11,OFF,6;
imgMask = []; % image Mask, [xrange; yrange];
wantdemean = 1; % whether to demean the data 
K = 2; % how many component you want

%%
%posiTuningImg = imgON(3:end); % Remove the condition -1 and 0
posiTuningImg = imgON(1:end); % Remove the condition -1 and 0
posiTuningImg = cellfun(@(x) permute(x,[1 2 4 3]), posiTuningImg, 'UniformOutput',0); % height X width X trial X time;
posiTuningImgCopy = posiTuningImg;

%% Preprocess the data for all condition
posiTuningImg = cellfun(@(x) squeeze(nanmean(x, 3)), posiTuningImgCopy, 'UniformOutput',0); % average trials
posiTuningImg = cellfun(@(x) x(:, :, timeWindow), posiTuningImg, 'UniformOutput',0); % extract time window;

if ~isempty(imgMask) % mask out pixels;
    for i=1:numel(posiTuningImg), posiTuningImg{i}(imgMask(1,1):imgMask(1,2), imgMask(2,1):imgMask(2,2), :)=nan; end %
end

posiTuningImg = cellfun(@(x) reshape(x,[316*316 length(timeWindow)]), posiTuningImg, 'UniformOutput',0); % reshape the image
% Exclude NaN pixels and transpose
nanMask = cellfun(@(x) ~any(isnan(x),2), posiTuningImg, 'UniformOutput',0);
posiTuningImg = cellfun(@(x,y) x(y,:)', posiTuningImg, nanMask, 'UniformOutput',0); % time X pixel*pixel in each cell
if wantdemean
    posiTuningImg = cellfun(@(x) x-repmat(mean(x,2),[1, size(x,2)]), posiTuningImg, 'UniformOutput',0);
end

%% --------------- do grand PCA --------------------- 
% We perform a PCA on the combined images from the average-trial images
% across all conditions

% combine all avg images
grandImg = cat(2,posiTuningImg{:});

% Run singular vector decomposition
% the U is the eigen vectors sorted by its eigen values
[U,S,VT] = svd(grandImg,'econ');

% let's visualize the first K components
close all;
h=figure;
set(h,'Position',[0 0 500 400]);
lh=myplot((timeWindow-stimOnset) * 0.2, -U(:,1:K)');
straightline(0,'v','k'); % add the stimulus onset line
xlabel('SOA (sec)');ylabel('Signal change (%)');
set(lh(1),'Color','r');
set(lh(2),'Color','b');
legend(lh,{'PC1','PC2'}, 'Box','off');
set(gca,'Color','none');


savefig(h,'testB_main_ON_PC.fig');

%% calculate the projected the image average trials for the 1st pc
compoGrand = -U(:,1); % We only project the first component
% We have to adjust the polarity of the PCs
valueLoad = cellfun(@(x) compoGrand'*x, posiTuningImg, 'UniformOutput',0); % pixel*pixel X 2
valueImg = repmat({NaN(1,316*316)},[1,numel(posiTuningImg)]);
valueImg = cellfun(@(x) reshape(x, [316 316]), valueImg, 'UniformOutput',0);
for i=1:numel(posiTuningImg); valueImg{i}(nanMask{i})=valueLoad{i}; end

valueImg_mat = cat(3, valueImg{:}); % concatenate images

% substract the mean images
%valueImg_mat_demean = valueImg_mat-nanmean(valueImg_mat,3);

valueImg_mat_demean = valueImg_mat;
h2 = figure;
imagesc(makeimagestack(valueImg_mat_demean)); caxis([-0.5, 0.5]);
colorbar();
saveas(h2,'~/Dropbox/stonesync/19attentionprobV1optimaging/B_mainExpON_PCA_load.png');