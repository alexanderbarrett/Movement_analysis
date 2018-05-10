% RANDOM FOREST SIGNATURE
function [candidateframes,Mdl_test] = findsimilarframes_mlfeatures_multifeatures_3(ML_features,labels,framesubset_search)


  
[testval,~] = Mdl_test.predict(ML_features(framesubset_search,:));
candidateframes = cellfun(@str2num,testval);

end