function [cv_error_output,confusion_matrix,cfmat] = cv_error(ML_features,cp,frames,labels,cp_num,fieldnames_beh)
    
    % create training and testing sets
    training_frames = frames(training(cp,cp_num));
    training_labels = labels(training(cp,cp_num));
    testing_frames = frames(test(cp,cp_num));

    % train model on training set, create prediction for testing set
    % RANDOM FOREST
%     [treebagger_labels,~,~] = ...
    % KNN
    [treebagger_labels,~] = ...
    findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);

    % calculate cross validation error
    num_correct = 0;
    hand_labels = labels(test(cp,cp_num));
    for n = 1:length(testing_frames)
       if hand_labels(n) == treebagger_labels(n)
           num_correct = num_correct + 1;
       end
    end
    cv_error_output = 1 - (num_correct / length(testing_frames));
    
    confusion_matrix = confusionmat(hand_labels,treebagger_labels);
    
    cfmat = zeros(length(fieldnames_beh) - 1);
    for n = 1:(length(fieldnames_beh) - 1)
        for i = 1:(length(hand_labels))
            if hand_labels(i) == n
                num = treebagger_labels(i);
                cfmat(n,num) = cfmat(n,num) + 1;
            end
        end
    end

end

