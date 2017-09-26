function [cluster_struct] = clusterMetricsTodd_JDM(cluster_struct)
%%

% dependences: https://github.com/de-Bivort-Lab/behavior-mapping/tree/master/metrics
% wtAll is (normalized!) wavelet data.
addpath('C:\Users\Jesse Marshall\Documents\GitHub\Movement_analysis\Common\metrics');



%%
numClusters = length(unique(cluster_struct.labels));

lb = 20; % # frames of short count
[xCount, shortCount, meanDwell, runLengths] = metricStateTransitions(cluster_struct.labels); %this input doesnt work

markovLLR = metricMarkovLLRatio(numClusters,cluster_struct.labels);

entropy = metricEntropy(cluster_struct.labels); % is this really informative?

exits = metricExitStates(numClusters,cluster_struct.labels);


    meanClusterDist = metricMeanClusterDist(numClusters,reshape(cluster_struct.labels,[],1),cluster_struct.wtAll);



cluster_struct.shortCount = shortCount;
cluster_struct.markovLLR = markovLLR;
cluster_struct.entropy = entropy;
cluster_struct.exits = exits;
cluster_struct.meanDwell = meanDwell;
cluster_struct.meanClusterDist = meanClusterDist;

end