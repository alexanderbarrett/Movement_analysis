% RANDOM FOREST SIGNATURE
function [candidateframes,score,Mdl_test,heldout_instances] = findsimilarframes_mlfeatures_multifeatures(ML_features,labels,framesubset_search,annotation_window)

% KNN SIGNATURE
% function [candidateframes,score] = findsimilarframes_mlfeatures_multifeatures(ML_features,labels,framesubset_search,annotation_window)


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


% annotation_window is frames_to_annotate
% featurepredict = featureuse(annotation_window,:);

                            use_knn = 0;
                            use_rf = 1;
                            
                            if (use_rf)
%   prior_vec = [0.06139; 0.00028; 0.03134; 0.01165; 0.02431; 0.03116; 0.00035; ...
%       0.00381; 0.00374; 0.00956; 0.01358; 0.02103; 0.01458; 0.0198; 0.06631; ...
%       0.00641; 0.03433; 0.00012; 0.02697; 0.18472; 0.03733; 0.16564; 0.23112; 0.00047];
%   
  %% collapse each instance of a behavior to an individual frame
  unique_behaviors = unique(labels);
  instances_labelled = [];
  features_indiv_instances = [];
  for beh_here = unique_behaviors
     behavioral_frames = find(labels == beh_here);
     aa = zeros(1,max(behavioral_frames));
     aa(behavioral_frames) = 1;
     labelledinstances = bwconncomp(aa);
     for mm = 1:numel(labelledinstances.PixelIdxList)
         %% include every 20th frame of a large instance, or take the average of a short instance
         numexamples = floor(numel(labelledinstances.PixelIdxList{mm})./20);
         if (numexamples <= 1)
         instances_labelled = cat(1,instances_labelled,beh_here);
         features_indiv_instances =cat(1,features_indiv_instances,...
             nanmean(featureuse_train(labelledinstances.PixelIdxList{mm},:),1));
         else
              instances_labelled = cat(1,instances_labelled,beh_here*ones(numexamples,1));
         features_indiv_instances =cat(1,features_indiv_instances,...
             (featureuse_train(randsample(labelledinstances.PixelIdxList{mm},numexamples),:)));
         end
     end
  end
  
  %% BLOCKED SAMPLE FOR CROSS VALIDATION
  holdout_frac = 0.9;
  num = floor((1 - holdout_frac) * length(instances_labelled)) ; 
 
    Mdl_test = TreeBagger(40,features_indiv_instances((1:num),:),...
        instances_labelled(1:num),'Method','classification','OOBPrediction','on');

% original_labels =     instances_labelled(heldoutsubset);

  
%   Mdl_test = TreeBagger(40,featureuse_train,labels,'Method','classification','OOBPrediction','on'); % , ...
%   Mdl_test = TreeBagger(40,features_indiv_instances,instances_labelled,'Method','classification','OOBPrediction','on');
  
                                                % 'MaxNumSplits',1000,'Prior',prior_vec);
 %   Mdl_test = TreeBagger(100,featureuse_train,labels,'Method','regression','OOBPrediction','on');

                           Mdl_test.oobError();
%     [testval,score,stddevs] = Mdl_test.predict(featurepredict);
    [testval,score,stddevs] = Mdl_test.predict(features_indiv_instances((num+1:end),:));
 candidateframes = cellfun(@str2num,testval);
 
%  cfmat = confusionmat(instances_labelled(heldoutsubset),candidateframes) ; 

heldout_instances = instances_labelled(num+1:end) ;
 
 
                            elseif (use_knn)
CVKNNMdl2 = fitcknn(featureuse_train,labels,'Distance','Euclidean',...
    'NumNeighbors',50,'Standardize',1);
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