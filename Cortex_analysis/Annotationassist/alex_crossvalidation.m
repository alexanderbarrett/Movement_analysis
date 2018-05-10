%%%%%
% 1. NEED TO CREATE A PARTITION :
%    c = cvpartition(group,'KFold',k) creates a random partition for a 
%    stratified k-fold cross validation. 
%    group is a numeric vector, categorical array, character array, or cell
%    array of character vectors indicating the class of each observation. 
%    Each subsample has roughly equal size and roughly the same class 
%    proportions as in group. 
%    cvpartition treats NaNs or empty character vectors in group as missing 
%    values.
%    k is the number of disjoint subsamples
%
% 2. FIND LOSS ERROR ESTIMATE
%    use crossval
%    possibly : mcr = crossval('mcr',X,y,'Predfun',predfun)
%%%%%

% CREATE THE PARTITION
% outputvector has labels every frame with a number, each number
% corresponds to a behavior. The behaviors are in fieldnames_beh. 
% see new_annotation_master_script.m line 34
group = outputvector;
c = cvpartition(group,'Kfold',10);

% FIND LOSS ERROR ESTIMATE
y = fieldnames_beh; 
classf= @(XTRAIN,ytrain,XTEST) nominal(classRF_predict(XTEST,classRF_train(XTRAIN,ytrain,1000)));
cvMCR = crossval('mcr',X,y,'predfun',classf,'partition',c)

y = candidate_frames;
X = outputvector';
cp = cvpartition(y,'k',10); % Stratified cross-validation

classf = @(XTRAIN, ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain));

cvMCR = crossval('mcr',X,y,'predfun',classf,'partition',cp)
