function compare_spectrogram_features(mocapstruct_pre,timeperiod_pre,mocapstruct_post,timeperiod_post)

figure(337)
marker_view = 1;

params.fps = mocapstruct.fps;
markers_pre_clipped = hipass_clip(mocapstruct_pre.markers_aligned_preproc,unique(cat(2,mocapstruct_pre.bad_frames_agg{marker_view},mocapstruct_pre.bad_frames_agg{5})),params);
markers_post_clipped = hipass_clip(mocapstruct_post.markers_aligned_preproc,unique(cat(2,mocapstruct_post.bad_frames_agg{marker_view},mocapstruct_post.bad_frames_agg{5})),params);

concat_trace = (cat(1,markers_pre_clipped.(mocapstruct_pre.markernames{marker_view})(:,:),...
    markers_post_clipped.(mocapstruct_pre.markernames{marker_view})(:,:)));

time_1 = 1:size(markers_pre_clipped.(mocapstruct_pre.markernames{marker_view}),1);
time_2 = bsxfun(@plus,1:size(markers_post_clipped.(mocapstruct_pre.markernames{marker_view}),1),(time_1(end)));



figure(333)
[~,fr_temp,time_clustering,pc_spectrograms_temp] = spectrogram(concat_trace(:,3),mocapstruct_pre.fps,...
     mocapstruct_pre.fps./2,1:1:50,mocapstruct_pre.fps);
 imagesc(log10(pc_spectrograms_temp))
 
     
        [coeff2,score2,latent2,tsquared2,explained2] = pca(log10(pc_spectrograms_temp(:,time_2_spect))','Centered','on');
        pca_agg_spectrogram =(score2');


        figure(2345)
        imagesc(coeff2)
        
        score2_zscored = bsxfun(@rdivide,bsxfun(@minus,score2,mean(score2,1)),std(score2,[],1));
        
        
        
time_1_spect = 1:floor(size(score2,1)*max(time_1)./max(time_2));
time_2_spect = (max(time_1_spect)+1):size(score2,1);        

        figure(33333)
        imagesc(score2_zscored')
        
   zscore_ratio =      nanmean(score2_zscored(time_1_spect,:),1)./...
            (abs(nanmean(score2_zscored(time_2_spect,:),1)));
        
        [vals,inds] = sort(abs(zscore_ratio),'Descend');
        
         figure(2346)
         subplot(2,1,1)
        imagesc(coeff2(:,inds))
        subplot(2,1,2)
       plot(vals)
        
        
        
end