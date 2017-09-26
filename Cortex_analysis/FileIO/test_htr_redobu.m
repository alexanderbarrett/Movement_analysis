
filearray = {'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid1.htr',...
    'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid2.htr',...
    'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid3.htr',...
    'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid4.htr',...
    'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid5.htr',...
    };


%'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid6.htr'

[limbposition,rotations,translations_local,translations_global] = concat_htr(filearray);

limbstart = 1;
limbnames = fieldnames(limbposition);
limbstop = numel(rotations );

%% visualize skeleton quickly
figure(23)
for mm = 1:30:10000
    for kk = limbstart:limbstop
   plot3( [squeeze(limbposition.(limbnames{kk})(1,1,mm)) ...
                        squeeze(limbposition.(limbnames{kk})(1,2,mm)) ],...
                    [squeeze(limbposition.(limbnames{kk})(2,1,mm)) ...
                        squeeze(limbposition.(limbnames{kk})(2,2,mm)) ],...
                       [squeeze(limbposition.(limbnames{kk})(3,1,mm)) ...
                        squeeze(limbposition.(limbnames{kk})(3,2,mm)) ],'linewidth',3)
                    hold on
    end
     
    hold off
    
       zlim([0 300])
            xlim([-400 400])
            ylim([-400 400])
                        set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
    pause(0.05)
end

rotations_eul = cell(1,numel(rotations));
for mm = limbstart:limbstop
rotations_eul{mm} = rotm2eul(rotations{mm});
end

%% bad frame detector -- this could be very slow, may just want to run for all at once
 bad_threshold = 5;
 agg_features = cat(1,diff(squeeze(cell2mat(translations_local')),1,2), diff(cell2mat(rotations_eul)',1,2));
 bad_frames = find(abs(agg_features)>bad_threshold);
 [~,bad_inds] = ind2sub(size(agg_features),bad_frames);
 bad_inds = unique(bad_inds);
 agg_features(:,bad_inds) = 0;
 agg_features(isnan(agg_features)) = 0;

  agg_features = bsxfun(@rdivide,bsxfun(@minus,agg_features,nanmean(agg_features,2)),nanstd(agg_features,[],2));
feat_mean = nanmean(agg_features,2);
goodind = find(~isnan(feat_mean));

 figure(23)
 imagesc(agg_features(goodind,:))
 caxis([-5 5])
 
 agg_features = agg_features(goodind,:);
% agg_features_global = cell2mat(transvelocity_global');
% bad_frames_global = find(abs(agg_features_global)>bad_threshold);
% agg_features_global(bad_frames_global) = 0;



%% do a quick clustering
%labels = WaveletCluster(agg_features',1:size(agg_features,2))
fps = 300;
num_pcs_1 = 20;
num_pcs_2 = 50;
num_clusters = 100;
clustering_ind = 1:size(agg_features,2);
clustering_overlap = fps./2;
clustering_window = fps;


[clusterobj,featurespectrogram] = Cluster_GMM(agg_features,fps,fps./2,fps,num_pcs_1,num_pcs_2,num_clusters);
%[clusterobj_global,featurespectrogram_global] = Cluster_GMM(agg_features_global,fps,fps./2,fps,num_pcs_1,num_pcs_2,num_clusters);

[~,sort_ja] = sort(clusterobj);
%[~,sort_global] = sort(clusterobj_global);

cluster_numbers = zeros(1,num_clusters);
for mm = 1:num_clusters
    cluster_numbers(mm) = numel(find(clusterobj == mm));
end

%% look at the features
figure(23)
imagesc(featurespectrogram(:,sort_ja))
caxis([-0.1 0.1])
% 
% figure(24)
% imagesc(featurespectrogram_global(:,sort_global))
% caxis([-2 2])

%% animate the skeleton
M = animate_skeleton(limbposition,1:10:1000,limbnames,limbstart,1)

               time_ordering_fulltrace = cell(1,num_clusters);

        for ll = 1:num_clusters

             [~,~,times_here,~] = spectrogram(squeeze(agg_features(1,:)),clustering_window,clustering_overlap,fps,fps);

            time_ordering_fulltrace{ll} = (unique(bsxfun(@plus,round(times_here(find(clusterobj==ll))*clustering_window)',...
                -floor(clustering_overlap):floor(clustering_overlap))));
            
            time_ordering_fulltrace{ll}(time_ordering_fulltrace{ll}<1) = 1;
            
            time_ordering_fulltrace{ll}(time_ordering_fulltrace{ll} > max( clustering_ind)) = max( clustering_ind);
          time_ordering_fulltrace{ll}=   clustering_ind(time_ordering_fulltrace{ll});
        end
        
        
        figure(37)
subplot(2,1,1)
test = cat(2,time_ordering_fulltrace{:});
imagesc(agg_features(:,test))
caxis([-0.1 0.1])

subplot(2,1,2)
plot(clusterobj(sort_ja))

        
        
frames_use = 10000;
savedirectory_subcluster = 'E:\Bence\Data\MOCAP\regression_animations\skeletontest2\';
mkdir(savedirectory_subcluster)
for jjj = find(cluster_numbers>5);%58;%36;%good_clusters(1:end)
    matlab_fr = 10;
    figure(370)
    frame_inds = (time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj}))))';
    
    M = [];
    movie_output_temp = animate_skeleton(limbposition,frame_inds',limbnames,limbstart,jjj)
    %animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
    
    v = VideoWriter(strcat(savedirectory_subcluster,'movie',num2str(jjj),'.mp4'),'MPEG-4');
    open(v)
    
    writeVideo(v,movie_output_temp)
    close(v)
end


%% visualize specific features
%  [~,~,times_here,spectrogram_here] = spectrogram(squeeze(delta_markers(3,6,clustering_ind,4)),clustering_window,clustering_overlap,fps,fps);
%         [~,~,~,spectrogram_here_2] = spectrogram(squeeze(delta_markers(3,8,clustering_ind,4)),clustering_window,clustering_overlap,fps,fps);
%         [~,~,~,spectrogram_here_3] = spectrogram(squeeze(delta_markers(5,6,clustering_ind,4)),clustering_window,clustering_overlap,fps,fps);
%           [~,~,~,spectrogram_here_4] = spectrogram(squeeze(markers_preproc.(marker_names{1})(clustering_ind,3)),clustering_window,clustering_overlap,fps,fps);
%    % xlimits = [1 500];
%             xlimits = [1 numel(cluster_vals)];
            uper_freq = 75;
figure(23)
            subplot(4,2,1)
            imagesc(log10(spectrogram_here(:,time_ordering)))
            title('head to spine')
            caxis([-8 1])
            ylim([0 50])
            xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,3)
            imagesc(log10(spectrogram_here_2(:,time_ordering)))
            title('hip to spine')
            caxis([-8 0])
            ylim([0 50])
            ylabel('Frequency [Hz]')
            xlabel('time ')
            xlim(xlimits)
            
            subplot(4,2,5)
            imagesc(log10(spectrogram_here_4(:,time_ordering)))
            title('head height')
            caxis([-8 1])
            ylim([0 50])
            xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,7)
            area(cluster_vals')
            %xlim([1 numel(cluster_vals)])
            ylabel('Cluster Number')
            xlim(xlimits)

