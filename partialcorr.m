% performing the partial correlation analysis

clear all;close all;clc

% setup
deviantProbOrder = [9 3 6 8 2 5 7 1 4];
deviantProb = [0.3478 0.0425 0.0425 0.0425 0.0425 0.0425 0.0425 0.0425 0.3478];

%% ========= ON exp ===============

% load the data from both monkeys
attmodu_B = load('mainExp_attmodu_B_ON.mat');
attmodu_K = load('mainExp_attmodu_K_ON.mat');


residual_hitRate_B = attmodu_B.hitRate - attmodu_B.hitRate/deviantProb * deviantProb;
residual_attModu_B = attmodu_B.betaAttInOutDiff(1:end-1) - attmodu_B.betaAttInOutDiff(1:end-1) / deviantProb  * deviantProb;

residual_hitRate_K = attmodu_K.hitRate - attmodu_K.hitRate/deviantProb * deviantProb;
residual_attModu_K = attmodu_K.betaAttInOutDiff(1:end-1) - attmodu_K.betaAttInOutDiff(1:end-1) / deviantProb  * deviantProb;

close all;
h1 = figure;
scatter(residual_hitRate_B, residual_attModu_B, 'r'); hold on;
scatter(residual_hitRate_K, residual_attModu_K, 'k'); hold on;
mycorrelation([residual_hitRate_B residual_hitRate_K], [residual_attModu_B, residual_attModu_K]);
xlabel('Residual hit rate');ylabel('Residual positional modulation (%)');

%close all; % plot the correlation of raw hitRate and attmodu
h2 = figure;
%scatter(residual_hitRate_B, residual_attModu_B, 'r'); hold on;
%scatter(residual_hitRate_K, residual_attModu_K, 'k'); hold on;
h=mycorrelation([attmodu_B.hitRate attmodu_K.hitRate], [attmodu_B.betaAttInOutDiff(1:end-1), attmodu_K.betaAttInOutDiff(1:end-1)]);
xlabel('Raw hit rate');ylabel('Raw positional modulation (%)');




%% ========= OFF exp ===============