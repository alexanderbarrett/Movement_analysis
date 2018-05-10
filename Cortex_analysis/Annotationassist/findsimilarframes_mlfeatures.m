function [candidateframes,score] = findsimilarframes_mlfeatures(ML_features,framelist,framesubset_search,annotation_window)


featureuse = cat(2,ML_features.pose_score(:,:),ML_features.trunk_vel(1,:)',ML_features.RGroom(:,1));
featureuse_train =featureuse(framesubset_search,:);

labels = zeros(1,numel(framesubset_search));
[~,ind_true_subset] = intersect(framesubset_search,framelist);
labels(ind_true_subset) = 1;
                            featurepredict = featureuse(annotation_window,:);

                            use_knn = 0;
                            use_rf = 1;
                            
                            if (use_rf)
  Mdl_test = TreeBagger(50,featureuse_train,labels,'Method','classification','OOBPrediction','on');
 %   Mdl_test = TreeBagger(100,featureuse_train,labels,'Method','regression','OOBPrediction','on');

                           Mdl_test.oobError();
    [testval,score,stddevs] = Mdl_test.predict(featurepredict);
 candidateframes = str2num(cell2mat(testval));
                            elseif (use_knn)
CVKNNMdl2 = fitcknn(featureuse_train,labels,'Distance','Euclidean',...
    'NumNeighbors',10,'Standardize',1);
[candidateframes,score] = predict(CVKNNMdl2,featurepredict);
                            end
  
%                             Mdl_test.error(regression_input(kfold_training_data([1:100000 120000:250000]) ,:),...
%                                 squeeze(regression_output(kfold_training_data([1:100000 120000:250000]),ll )));
%                             
%                             Mdl_test.error(regression_input(kfold_training_data(100000:110000) ,:),...
%                                 squeeze(regression_output(kfold_training_data(100000:110000),ll )));
%                             
%                             testval = Mdl_test.predict(regression_input(kfold_training_data(1:100000) ,:));
%                             sqrt(mean( ( squeeze(regression_output(kfold_training_data(1:100000),ll ))-testval).^2))
%                             
%                            
%                             sqrt(mean( ( squeeze(regression_output(kfold_training_data(105000:115000),ll ))-testval).^2))


end