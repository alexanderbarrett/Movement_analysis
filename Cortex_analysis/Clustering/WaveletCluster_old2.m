function [labels , feature_mu, feature_sigma,wtAll,fr]= WaveletCluster_old(features, frames, opts)

% DEPENDENCIES:
% addpath(genpath(...\Movement_analysis-master\MotionMapper-master\));
% addpath(genpath(...\Movement_analysis-master\VbGm));

% clustering script using wavelet space, and GMM / vb method
%
% features is T x d (observations x dimensions)
%  - this assumes d is relatively small (~<50)!!!
%

%% initialize options
if nargin<3;
    opts.whiten = 0;
    opts.frameNormalize = 0;
    opts.clustermethod = 'GMM';
    num = size(features,1); % number modes (spectrograms to find) (usually want full dimension)
    % num = pcuse;
    ds = 1; % down sampling
	samprate = 100;
    params = struct;
    params.samplingFreq = samprate/ds;
    params.numPeriods=25; %distinct number of frequencies to use
    params.minF = 0.6; % min freq to analyze
    params.maxF = 50; % max freq to analyze
    
    
    % for gmm model parameters
    pcuse = 20;
    numclusters = 100;
    lambda = 0.1; % regularization
else
    num=opts.num;
    ds = opts.ds;
    samprate = opts.samprate;
    params = opts.params;
    
    pcuse = opts.pcuse;
    numclusters = opts.numclusters;
    lambda = opts.lambda;
end

    
    subsample = 10; % number of times to run GMM fit on subsample
    numsample = 500; % number of data points to sample
    


%% crop frames
if ~isempty(frames);
features = features(frames,:);
end
T = size(features,1);
d = size(features,2);

%% whiten the data
if opts.whiten
features = bsxfun(@rdivide,bsxfun(@minus,features,mean(features,2)),std(features,[],2));
end

%% find wavelet transform (dependent on berman code!!!)
disp('finding wavelet transform');
[wt,fr]=findWavelets(features,num,setRunParameters(params));

% wt should be T x d*numPeriods

if opts.frameNormalize;

% frame normalize (see Todd paper)
wtAmps = sum(wt,2);
% exclude low energy frames
options = statset('MaxIter',1000);
gmmodel = fitgmdist(wtAmps,5,'Options',options,'replicates',1);
% thresh as min between two largest peaks
range = quantile(wtAmps,[.01 .99]);
range = linspace(range(1)-std(wtAmps),range(2)+std(wtAmps),1000);
gaus_fit = pdf(gmmodel,range');
[pks,loc] = findpeaks(gaus_fit);
[~,threshind] = min(gaus_fit(loc(1):loc(2))); threshind = threshind + loc(1);
% figure;plot(range,gaus_fit);
% hold on;plot(range(loc),pks,'ro');
% plot(range(threshind),gaus_fit(threshind),'go');

thresh = range(threshind);
useind = find(wtAmps > thresh);
keepind = find(wtAmps < thresh);

%figure; hold on;
%histogram(wtAmps(wtAmps > thresh)); histogram(wtAmps(wtAmps < thresh));

% frame normalize
wtNorm = bsxfun(@rdivide,wt,wtAmps);
wtNorm = wtNorm .* max(max(wt))/max(max(wtNorm)); % rescale
% figure;imagesc(wtNorm(:,1:100)');colorbar; title('wt norm');
% figure;imagesc(wt(:,1:100)');colorbar; title('wt raw');
wtAll = wtNorm;
wtAll(keepind,:) = wt(keepind,:);
% figure;imagesc(wtAll(:,1:100)');colorbar; title('wt all');
else
wtAll = wt; % because frame normalization still a bit sketchy
end


%% take PCA of wavelets
disp('finding PCs');

[U,S,V] = pca_randomized(wtAll, pcuse);
%[coeff,score,latent,tsquared,explained,mu] = pca(wtAll);
coeff = zeros(size(wtAll,2),size(wtAll,2));

score = 1000*U;
coeff = S;

%% cluster PC space
disp('clustering');




if strcmp(opts.clustermethod,'GMM');

% gmmodel
options = statset('MaxIter',1000,'Display','final');
% subsample = 10;
% numsample = 500;
fit = {}; nlogl_all=[];
for q = 1:subsample
    subInd = randsample(size(score,1),numsample);
    %fit = fitgmdist(score(:,1:pcuse), numclusters,'Options',options,'RegularizationValue',lambda);
    %[labels,nlogl,P,logpdf,M] = cluster(fit,score(:,1:pcuse));
    fit{q} = fitgmdist(score(subInd,1:pcuse), numclusters,'Options',options,'RegularizationValue',lambda);
    [~,nlogl] = cluster(fit{q},score(:,1:pcuse));
    nlog_all(q) = nlogl;
end

[~,usefit] = min(nlog_all);
[labels,nlogl,P,logpdf,M] = cluster(fit{usefit},score(:,1:pcuse));

elseif strcmp(opts.clustermethod,'VB');

% vb ( see bishop page 474)
%addpath('C:\Users\Kevin\Documents\Movement_analysis\VbGm');
prior.alpha = .1; % Dirichlet: 'effective prior number of observations associated with each component of the mixture'
prior.kappa = 1; % Gaussian
prior.m = mean(score(:,1:pcuse)',2); % Gaussian (change to kmeans initialization?)
prior.v = d+1; %Wishart
prior.M = eye(d); %Wishart
[labels, model, L] = mixGaussVb(score(:,1:pcuse)',numclusters);
end

%% temp
[~,ind] = sort(labels);
%figure;imagesc(wtAllcrop(ind,:)');


%% convert fit back to wavelet dimension
% to account for full pca space

mu = fit{usefit}.mu;
sigma = fit{usefit}.Sigma;

feature_mu = mu * V';
feature_sigma = zeros(size(wtAll,2),size(wtAll,2) , numclusters);
for j = 1:numclusters
feature_sigma(:,:,j) = V * sigma(:,:,j) * V';
end

%% **** plotting
% plots wavelet feature space sorted into clusters
[sorted_labels,ind] = sort(labels,'ascend');
figure;imagesc(wtAll(labels,:));

label_ind = find(diff(sorted_labels)>0);
hold on;
for j = 1:length(label_ind);
    line([0 size(wtAll,2)],[label_ind(j) label_ind(j)],'color','r');
end

end