function predict_markers_position(mocapstruct,mod_1)



struct_fn = fieldnames(mocapstruct.markers_aligned_preproc);
clustering_type =2;
feature_mat_1 = [];

plotdir_here = strcat(mocapstruct.plotdirectory,filesep,'comparison_subcluster',num2str(clustering_type),filesep);
mkdir(plotdir_here);

downsample = 3;

for mm = mod_1.cluster_markersets{clustering_type}
    
    frame_ex_1 = intersect(mod_1.clustering_inds_agg{clustering_type}(1:downsample:end),mocapstruct.move_frames);
    
    feature_mat_1 = cat(2,feature_mat_1,mocapstruct.markers_aligned_preproc.(struct_fn{mm})(frame_ex_1,:));
end


%% compare a subgroup




do_subset = 0;
%% first choose which groups to compare
if (do_subset)
fit_ind = 1:500:1000000;
marker_groups = {[1:3],[4,6,7,8],[8,9,18],[8,10,17]};

comparison_groups = [1,2,3,4,5,7,8,9,18,10,17] ;
fitlosses = zeros(numel(comparison_groups),numel(comparison_groups));
fitlosses_norm = zeros(numel(comparison_groups),numel(comparison_groups));

for group_1 = comparison_groups
    for group_out = setxor(comparison_groups,group_1)
        fprintf('for group_1 %f group_2 %f \n',group_1,group_out);
        
        markers_in = find(group_1 == mod_1.cluster_markersets{clustering_type} );
        markers_out = find(group_out == mod_1.cluster_markersets{clustering_type} );
        
        index_in = 3*(markers_in-1)+1:3*(markers_in);
        % loop over all markers in the output group
        
        
        fit_in = feature_mat_1(fit_ind,index_in);
        mdlloss = 0;
        mdlloss_norm = 0;
        
        for jj = 1:3
            fit_out = feature_mat_1(fit_ind,3*(markers_out-1)+1);
            
            Mdl = fitrlinear( fit_in,fit_out,'KFold',10);
            mdlloss =mdlloss+ Mdl.kfoldLoss;
            mdlloss_norm =mdlloss_norm+ Mdl.kfoldLoss./mean(bsxfun(@minus,fit_out,mean(fit_out,1)).^2);
            
        end
        
        
        
        
        fitlosses(find(comparison_groups == group_1),find(comparison_groups == group_out)) = mdlloss./3;
        fitlosses_norm(find(comparison_groups == group_1),find(comparison_groups == group_out)) = mdlloss_norm./3;
        
        %Mdl = fitrsvm( fit_in,fit_out,'KFold',10,'KernelFunction','gaussian')
        
    end
end

fitlosses_max = max(fitlosses);
fitlosses_norm2 = bsxfun(@rdivide,fitlosses,fitlosses_max);

figure(5)
imagesc(   fitlosses_norm)
set(gca,'XTickLabels',mocapstruct.markernames( comparison_groups),'YTickLabels',mocapstruct.markernames( comparison_groups));


else
%% do all intermarker comparisons
fit_ind = 1:500:1000000;

%% first choose which groups to compare
comparison_groups = [1:4,6:20] ;
fitlosses = zeros(numel(comparison_groups),numel(comparison_groups));
fitlosses_norm = zeros(numel(comparison_groups),numel(comparison_groups));

for group_1 = comparison_groups
    for group_out = setxor(comparison_groups,group_1)
        fprintf('for group_1 %f group_2 %f \n',group_1,group_out);
                
        downsample = 3;
        bad_frames_here = cat(2,mocapstruct.bad_frames_agg{group_1},mocapstruct.bad_frames_agg{group_out},...
            mocapstruct.bad_frames_agg{5},mocapstruct.bad_frames_agg{6});
        frame_ex_1 = intersect(setxor(1:size(mocapstruct.aligned_mean_position,1),bad_frames_here),mocapstruct.move_frames);
        frame_ex_1 = frame_ex_1(1:500:min(numel(frame_ex_1),1000000));
        fit_in = mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{group_1})(frame_ex_1,:);
        
        mdlloss = 0;
        mdlloss_norm = 0;
        for jj = 1:3
            fit_out = mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{group_out})(frame_ex_1,jj);
                        
            Mdl = fitrlinear( fit_in,fit_out,'KFold',10);
            mdlloss =mdlloss+ Mdl.kfoldLoss;
            mdlloss_norm =mdlloss_norm+ Mdl.kfoldLoss./mean(bsxfun(@minus,fit_out,mean(fit_out,1)).^2);            
        end

        fitlosses(find(comparison_groups == group_1),find(comparison_groups == group_out)) = mdlloss./3;
        fitlosses_norm(find(comparison_groups == group_1),find(comparison_groups == group_out)) = mdlloss_norm./3;
        
        %Mdl = fitrsvm( fit_in,fit_out,'KFold',10,'KernelFunction','gaussian')
        
    end
end

fitlosses_max = max(fitlosses);
fitlosses_norm2 = bsxfun(@rdivide,fitlosses,fitlosses_max);

figure(5)
imagesc(   fitlosses_norm)
caxis([0 1.2])
xlabel('Predicted Marker')
ylabel('Predictor Marker')
h=colorbar;
ylabel(h, 'Cross-Val MSE/MS Power')
set(gca,'XTick',1:numel(comparison_groups),'YTick',1:numel(comparison_groups),'XTickLabels',mocapstruct.markernames( comparison_groups),'YTickLabels',mocapstruct.markernames( comparison_groups));
title('Marker modularity')
print('-dpng',strcat(mocapstruct.plotdirectory,filesep,'markermodularity.png'))


Z = linkage(fitlosses_norm,'single');
figure(11)
dendrogram(Z,'labels',mocapstruct.markernames( comparison_groups))
ylabel('Nearest Neighbor distance between clusters')
print('-dpng',strcat(mocapstruct.plotdirectory,filesep,'markerhierarchy.png'))



         figure(4)
    plot(fit_out(fitpts,jj),'r')
    hold on
    plot( Mdl.kfoldPredict,'g')
    hold off
%         

end


%% do dynamices as well

%% first choose which groups to compare
%comparison_groups = [1:4,6:20] ;
comparison_groups = [1:4,6:20] ;

fitlosses = zeros(numel(comparison_groups),numel(comparison_groups));
fitlosses_norm = zeros(numel(comparison_groups),numel(comparison_groups));

for group_1 = comparison_groups
    for group_out = setxor(comparison_groups,group_1)
        fprintf('for group_1 %f group_2 %f \n',group_1,group_out);
                
%% get the input spectrograms, clipped for each
        bad_frames_here = cat(2,mocapstruct.bad_frames_agg{group_1},mocapstruct.bad_frames_agg{group_out},...
            mocapstruct.bad_frames_agg{5},mocapstruct.bad_frames_agg{6});
        
        fragment_frames = intersect(setxor(1:size(mocapstruct.aligned_mean_position,1),bad_frames_here),mocapstruct.move_frames);
        params.fps = 300;
        markers_here = struct('singlemarker',[]);
     markers_here.singlemarker = mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{group_1});
         [ markers_here,clipped_index] =  hipass_clip_fragments(markers_here,fragment_frames,params);
       % fit_in = markers_here.singlemarker;
        
        
        
opts.fps =300./downsample;
opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;
      [dyad_out,fr,time_clustering] =   get_dyadic_spectrogram(markers_here.singlemarker',opts);
      
[COEFF, SCORE_dyn, LATENT, TSQUARED,EXPLAINED] = pca(dyad_out');
      fit_in = SCORE_dyn(:,1:10);

        
      %% do for output as well
       markers_here = struct('singlemarker',[]);
     markers_here.singlemarker = mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{group_out});
         [ markers_here,clipped_index] =  hipass_clip_fragments(markers_here,fragment_frames,params);
               [dyad_out2,fr,time_clustering] =   get_dyadic_spectrogram(markers_here.singlemarker',opts);
[COEFF, SCORE_dyn2, LATENT, TSQUARED,EXPLAINED] = pca(dyad_out2');
      fit_out = SCORE_dyn2(:,1:3);

         
      %  frame_ex_1 = 
        
        
       % frame_ex_1 = frame_ex_1(1:500:min(numel(frame_ex_1),1000000));
     %   fit_in = mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{group_1})(frame_ex_1,:);
        
     fitpts = randsample(1:size(fit_out,1),min(20000,size(fit_out,1)));
     
        mdlloss = 0;
        mdlloss_norm = 0;
        for jj = 1:3
          %  fit_out = mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{group_out})(frame_ex_1,jj);
                        
            Mdl = fitrlinear( fit_in(fitpts,:),fit_out(fitpts,jj),'KFold',10);
            mdlloss =mdlloss+ Mdl.kfoldLoss;
            mdlloss_norm =mdlloss_norm+ Mdl.kfoldLoss./mean(bsxfun(@minus,fit_out(fitpts,jj),mean(fit_out(fitpts,jj),1)).^2);            
        end

        fitlosses(find(comparison_groups == group_1),find(comparison_groups == group_out)) = mdlloss./3;
        fitlosses_norm(find(comparison_groups == group_1),find(comparison_groups == group_out)) = mdlloss_norm./3;
        
        %Mdl = fitrsvm( fit_in,fit_out,'KFold',10,'KernelFunction','gaussian')
        
    end
end

fitlosses_max = max(fitlosses);
fitlosses_norm2 = bsxfun(@rdivide,fitlosses,fitlosses_max);

figure(50)
imagesc(   fitlosses)
%caxis([0.3 1])
xlabel('Predicted Marker')
ylabel('Predictor Marker')
h=colorbar;
ylabel(h, 'Cross-Val MSE/MS Power')
set(gca,'XTick',1:numel(comparison_groups),'YTick',1:numel(comparison_groups),'XTickLabels',mocapstruct.markernames( comparison_groups),'YTickLabels',mocapstruct.markernames( comparison_groups));
title('Marker modularity')
print('-dpng',strcat(mocapstruct.plotdirectory,filesep,'markermodularity_dyn.png'))


Z = linkage(fitlosses_norm,'single');
figure(110)
dendrogram(Z,'labels',mocapstruct.markernames( comparison_groups))
ylabel('Nearest Neighbor distance between clusters')
print('-dpng',strcat(mocapstruct.plotdirectory,filesep,'markerhierarchy_dyn.png'))









%

end