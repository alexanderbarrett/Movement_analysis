function [candidateframes,score,Mdl_test] = findsimilarframes_mlfeatures_multifeatures(ML_features,labels,framesubset_search,annotation_window)

% 
% featureuse = cat(2,ML_features.pose_score(:,:),ML_features.trunk_vel(1,:)',ML_features.head_vel(1,:)',ML_features.RGroom(:,1),...
%    ML_features.absolute_velocity_trunk_abs_100,ML_features.absolute_velocity_trunk_abs_100,...
%    ML_features.rel_velocity_hipL_abs_100,ML_features.rel_velocity_hipR_abs_100,...
%    ML_features.rel_velocity_head_abs_100,ML_features.rel_std_velocity_head_abs_100);
% 

featureuse=[];
%% load in features
[ tsne_features,num_feat]= load_tsne_features('annotation_ii');
    tsnefeat_name = cell(0,1);
   
        for mm = 1:numel(tsne_features)
              if mm>numel(num_feat)
          num_feat(mm) = 1;
              end
      feature_temp = ML_features.(tsne_features{mm});
      
      %% issue with transposition
      if size(feature_temp,1)<size(feature_temp,2)
         feature_temp = feature_temp'; 
      end
      
      %% issue with wavelet sizing
      if size(feature_temp,1)<numel(ML_features.framelist_true)
      feature_temp = cat(1,feature_temp,zeros(numel(ML_features.framelist_true)-size(feature_temp,1),size(feature_temp,2)));
      end
          
           featureuse = cat(2,featureuse ,feature_temp(:,1:num_feat(mm)));
        end
         for jj = 1:num_feat(mm)
            tsnefeat_name{end+1} = strcat(tsne_features{mm},'__',num2str(jj));
            end
 
%% training and testing features
featureuse_train =featureuse(framesubset_search,:);
% 
% labels = zeros(1,numel(framesubset_search));
% [~,ind_true_subset] = intersect(framesubset_search,framelist);
% labels(ind_true_subset) = 1;
                            featurepredict = featureuse(annotation_window,:);

                            use_knn = 0;
                            use_rf = 1;
                            
                            if (use_rf)
  Mdl_test = TreeBagger(50,featureuse_train,labels,'Method','classification','OOBPrediction','on');
 %   Mdl_test = TreeBagger(100,featureuse_train,labels,'Method','regression','OOBPrediction','on');

                           Mdl_test.oobError();
    [testval,score,stddevs] = Mdl_test.predict(featurepredict);
 candidateframes = cellfun(@str2num,testval);
 
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