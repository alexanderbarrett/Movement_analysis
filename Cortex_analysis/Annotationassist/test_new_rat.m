function [candidateframes,instances_labelled] = test_new_rat(ML_features,labels,framesubset_search,Mdl_test)

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
  
  %% RANDOMLY SAMPLE FRAMES for cross validation
%   holdout_frac = 0.3;
%  random_subset = randsample(1:numel(instances_labelled),floor((1-holdout_frac)*numel(instances_labelled)));
%     Mdl_test = TreeBagger(40,features_indiv_instances(random_subset,:),...
%         instances_labelled(random_subset),'Method','classification','OOBPrediction','on');
%     
% heldoutsubset = setxor(1:numel(instances_labelled),random_subset);
% original_labels =     instances_labelled(heldoutsubset);

  
%   Mdl_test = TreeBagger(40,featureuse_train,labels,'Method','classification','OOBPrediction','on'); % , ...
%   Mdl_test = TreeBagger(40,features_indiv_instances,instances_labelled,'Method','classification','OOBPrediction','on');
  
                                                % 'MaxNumSplits',1000,'Prior',prior_vec);
 %   Mdl_test = TreeBagger(100,featureuse_train,labels,'Method','regression','OOBPrediction','on');

%                            Mdl_test.oobError();
%     [testval,score,stddevs] = Mdl_test.predict(featurepredict);

% [testval,score,stddevs] = Mdl_test.predict(features_indiv_instances(heldoutsubset,:));
[testval,~] = Mdl_test.predict(features_indiv_instances);
 candidateframes = cellfun(@str2num,testval);
 
%  cfmat = confusionmat(instances_labelled(heldoutsubset),candidateframes) ; 


 
 

end