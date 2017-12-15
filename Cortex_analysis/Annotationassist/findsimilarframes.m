function candidateframes = findsimilarframes(mocapstruct,framelist,markerlist,framesubset_search)

fieldnames_here = fieldnames(mocapstruct.markers_aligned_preproc);

params.fps = 300;


opts.fps =300./1;
opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;

%% get the frames to serve as a template


%% get total features to search through
   %% do a hipass clip+interpolation to correct
          feature_mat_1 = [];
                    feature_mat_template = [];

     clipped_pre =  hipass_clip_fragments(mocapstruct.markers_aligned_preproc,framesubset_search,params);
markerout = hipass_clip_fragments(mocapstruct.markers_aligned_preproc ,framelist,params);

         for mm = markerlist
     feature_mat_1 = cat(2,feature_mat_1,clipped_pre.(fieldnames_here{mm}));
          feature_mat_template = cat(2,feature_mat_template,markerout.(fieldnames_here{mm}));

         end
nummarkershere = size(feature_mat_1,2)./3;

agg_features = permute(reshape(feature_mat_1,size(feature_mat_1,1),3,[]),[1 3 2]);
agg_features_template = permute(reshape(feature_mat_template,size(feature_mat_template,1),3,[]),[1 3 2]);

[dyadic_spectrograms_1,fr,timespect] = get_dyadic_spectrogram(feature_mat_1(:,:)',opts);
[dyadic_spectrograms_template,fr,timespect_template] = get_dyadic_spectrogram(feature_mat_template(:,:)',opts);

%% whiten each spectrogram vector
dyadic_spectrograms_template_w = bsxfun(@rdivide,bsxfun(@minus,dyadic_spectrograms_template,mean(dyadic_spectrograms_1,2)),std(dyadic_spectrograms_1,[],2));
dyadic_spectrograms_1_w = bsxfun(@rdivide,bsxfun(@minus,dyadic_spectrograms_1,mean(dyadic_spectrograms_1,2)),std(dyadic_spectrograms_1,[],2));

    
%mean_template = nanmean(dyadic_spectrograms_template_w,2);

%template_score = dyadic_spectrograms_1_w'*mean_template;

%[~,most_inform] = sort(mean_template,'DESCEND');
%dim_use = most_inform(1:50);
%template_score_subset = dyadic_spectrograms_1_w(dim_use,:)'*mean_template(dim_use);


distance_template = 
distance_std 


CVKNNMdl2 = fitcknn(X,Y,'Distance','Euclidean',...
    'NumNeighbors',k,'KFold',10,'Standardize',1);


figure(55)
subplot(3,1,1)
plot(mean_template)
subplot(3,1,2)
plot(template_score)
subplot(3,1,3)
plot((template_score_subset))

thresh = input('enter threshold');
candidateframes = 300*timespect(find((cat(1,template_score_subset))>thresh));
candidateframes = unique(bsxfun(@plus,candidateframes',-150:150));
candidateframes(candidateframes<1) = 1;

end

