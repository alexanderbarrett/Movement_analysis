%%%%%%%
% THIS FILE CONTAINS SCRATCHWORK
% CV ERROR FILE IS ''alex_cv2.m'
%%%%%%%


load('fisheriris');
y = species;
X = meas;
cp = cvpartition(y,'k',10); % Stratified cross-validation

classf = @(XTRAIN, ytrain,XTEST)(classify(XTEST,XTRAIN,...
ytrain));

cvMCR = crossval('mcr',X,y,'predfun',classf,'partition',cp)




data = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];

cp = cvpartition(data,'kfold',3)

a = training(cp,2)
data(a)
data(test(cp,2))


fields = fieldnames(all_behaviors);
for i = 1:numel(fields)
  fprintf(fields{i})
  fprintf('\n')
  display(all_behaviors.(fields{i})) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLOCKED PARTITION - TRAIN ON FIRST 90%, TEST ON LAST 10%

frames = annotated_frames_goodtracking;
labels = annotation_output_goodtracking;

% train on first 90%, test on last 10%  
cutoff = floor(0.9 * length(frames));

% create training and testing sets
training_frames = frames(1:cutoff);
training_labels = labels(1:cutoff);
testing_frames = frames((cutoff + 1):end);

% train model on training set, create prediction for testing set
[treebagger_labels,~,~] = ...
findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);

% calculate cross validation error
num_correct = 0;
hand_labels = labels((cutoff + 1):end);
for n = 1:length(testing_frames)
   if hand_labels(n) == treebagger_labels(n)
       num_correct = num_correct + 1;
   end
end
cv_error = 1 - (num_correct / length(testing_frames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BLOCKED K-FOLD CV ERROR

num_folds = 10;

frames = annotated_frames_goodtracking;
labels = annotation_output_goodtracking;

% need to drop out up to num_folds number of frames off of end 
% because reshape must produce a valid matrix
num_to_drop = mod(length(frames),num_folds);
frames = frames(1:(end - num_to_drop));
labels = labels(1:(end - num_to_drop));

% create the blocked partition
frames_partition = reshape(frames,[],num_folds);
labels_partition = reshape(labels,[],num_folds);

% initialization
indices = 1:num_folds;
blocked_cv_errors = zeros(1,num_folds);

for i = 1:num_folds
    
    % determine which part to hold out
    training_items = indices(indices~=i);
    
    % initialization
    training_labels = [0];
    training_frames = [0];
    
    % create vector for training set
    for j = 1:(num_folds - 1)
        j2 = training_items(j);
        training_labels = [training_labels,labels_partition(:,j2)'];
        training_frames = [training_frames,frames_partition(:,j2)'];              
    end
    
    % create vector for testing set
    testing_frames = frames_partition(:,i);
    
    % run treebagger
    [treebagger_labels,~,~] = ...
    findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);
    
    % calculate cross validation error
    num_correct = 0;
    hand_labels = labels_partition(:,i);
    for n = 1:length(testing_frames)
       if hand_labels(n) == treebagger_labels(n)
           num_correct = num_correct + 1;
       end
    end
    blocked_cv_error = 1 - (num_correct / length(testing_frames));
    
    % add cv error for this run to array
    blocked_cv_errors(i) = blocked_cv_error;
    
end  

% calculate avg cv error
blocked_cv_error_avg = mean(blocked_cv_errors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% train on JDM25 Caff, test on Vicon8 caff
JDM25_frames = annotated_frames_goodtracking;
JDM25_hand_labels = annotation_output_goodtracking;

% create training and testing sets
training_frames = JDM25_frames;
training_labels = JDM25_hand_labels;
testing_frames = Vicon8_frames;

% train model on training set, create prediction for testing set
[treebagger_labels,~,~] = ...
findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);

% calculate cross validation error
num_correct = 0;

for n = 1:length(testing_frames)
   if Vicon8_hand_labels(n) == treebagger_labels(n)
       num_correct = num_correct + 1;
   end
end
cv_error_output = 1 - (num_correct / length(testing_frames));


















