function [shortCount, markovLLR, entropy, exits, meanDwell, meanClusterDist] = clusterMetricsTodd(labels, wtAll);
%%

% dependences: https://github.com/de-Bivort-Lab/behavior-mapping/tree/master/metrics
% wtAll is (normalized!) wavelet data.
addpath('C:\Users\Kevin\Documents\Lab\behavior-mapping\metrics');



%%
numClusters = length(unique(labels));

lb = 20; % # frames of short count
[xCount, shortCount, meanDwell, runLengths] = metricStateTransitions(labels, lb);

markovLLR = metricMarkovLLRatio(numClusters,labels);

entropy = metricEntropy(labels); % is this really informative?

exits = metricExitStates(numClusters,labels);

if nargin < 2
    meanClusterDist = [];
else;
    meanClusterDist = metricMeanClusterDist(numClusters,labels',wtAll);
end


end