function [output_struct] = compare_features(mocapstruct_pre,mod_1,mocapstruct_post,mod_2,comptype,framesin_1,framesin_2,clustering_type,plotflag)

    struct_fn = fieldnames(mocapstruct_pre.markers_aligned_preproc);
    %if clustering type isn't an input
    if nargin< 8
    clustering_type =2;
    end
    
    feature_mat_1 = [];
        feature_mat_2 = [];

        plotdir_here = strcat(mocapstruct_pre.plotdirectory,filesep,'comparison_subcluster',num2str(clustering_type),filesep);
        mkdir(plotdir_here);
        
    downsample = 3;

    if (strcmp(comptype,'pose'))

    for mm = mod_1.cluster_markersets{clustering_type}    
        
 frame_ex_1 = intersect(mod_1.clustering_inds_agg{clustering_type}(1:downsample:end),mocapstruct_pre.move_frames);
    frame_ex_2 = intersect(mod_2.clustering_inds_agg{clustering_type}(1:downsample:end),mocapstruct_post.move_frames);
     frames_total = cat(2,frame_ex_1,frame_ex_2);

     
     if nargin> 5
         frame_ex_1 = intersect(mod_1.clustering_inds_agg{clustering_type}(1:downsample:end),framesin_1);
                  frame_ex_2 = intersect(mod_2.clustering_inds_agg{clustering_type}(1:downsample:end),framesin_2);

     end
     
     feature_mat_1 = cat(2,feature_mat_1,mocapstruct_pre.markers_aligned_preproc.(struct_fn{mm})(frame_ex_1,:));
          feature_mat_2 = cat(2,feature_mat_2,mocapstruct_post.markers_aligned_preproc.(struct_fn{mm})(frame_ex_2,:));
    end
    
    else
        %% do a hipass clip+interpolation to correct
          feature_mat_1 = [];
        feature_mat_2 = [];
 frame_ex_1 = intersect(mod_1.clustering_inds_agg{clustering_type}(1:end),mocapstruct_pre.move_frames);
    frame_ex_2 = intersect(mod_2.clustering_inds_agg{clustering_type}(1:end),mocapstruct_post.move_frames);
     
     if nargin> 5
         frame_ex_1 = intersect(intersect(mod_1.clustering_inds_agg{clustering_type}(1:end),framesin_1),mocapstruct_pre.move_frames);
                  frame_ex_2 = intersect(intersect(mod_2.clustering_inds_agg{clustering_type}(1:end),framesin_2),mocapstruct_post.move_frames);

     end
    
     
     
     params.fps = 300;
     
     clipped_pre =  hipass_clip_fragments(mocapstruct_pre.markers_aligned_preproc,frame_ex_1,params,mod_1.cluster_markersets{clustering_type});
          clipped_post =  hipass_clip_fragments(mocapstruct_post.markers_aligned_preproc,frame_ex_2,params,mod_1.cluster_markersets{clustering_type});

         for mm = mod_1.cluster_markersets{clustering_type}    

     feature_mat_1 = cat(2,feature_mat_1,clipped_pre.(struct_fn{mm}));
     feature_mat_2 = cat(2,feature_mat_2,clipped_post.(struct_fn{mm}));
    end
        
        
        
    end
    
    
             % feature_mat_2 = feature_mat_1(250000:end,:);

    % feature_mat_1 = feature_mat_1(1:250000,:);
if (strcmp(comptype,'pose'))


%% first compute PCS of the joint and look at sample 2 and 3d point clouds
feature_time_ind = 1;
%pca on both
%[COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(cat(feature_time_ind,feature_mat_1,feature_mat_2));
%pca on 1 to get rid of pcs that explain subtle differences
[COEFF, SCORE_1, LATENT, TSQUARED,EXPLAINED] = pca(feature_mat_1);

SCORE_2 = bsxfun(@minus,feature_mat_2,mean(feature_mat_2,1))*COEFF;
SCORE = cat(1,SCORE_1,SCORE_2);

maxinds = min(size(COEFF,1),15);

%% plot the PCs
num_colors = numel(mocapstruct_post.markercolor);

if (plotflag)
    h=figure(665)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'w'},num_colors ,1),mocapstruct_post.links,h)
view([-14 28])

  h=figure(669)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'w'},num_colors ,1),mocapstruct_post.links,h)
hold on
plot_eigenpose_subset(reshape(nanmean(feature_mat_2(:,:),1),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'r'},num_colors ,1),mocapstruct_post.links,h)
hold off
ylim([-100 200])
view([-14 28])

    
    
    
    
h=figure(666)
for ind_pick = 1:4%maxinds
subplot(1,4,ind_pick)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'w'},num_colors ,1),mocapstruct_post.links,h)
hold on
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])'+50.*reshape(COEFF(:,ind_pick),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'r'},num_colors ,1),mocapstruct_post.links,h)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])'-50.*reshape(COEFF(:,ind_pick),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'b'},num_colors ,1),mocapstruct_post.links,h)
view([-33 28])
hold off
end
end



feature_1_inds = 1:size(feature_mat_1,feature_time_ind);
feature_2_inds = (max(feature_1_inds)+1):(max(feature_1_inds)+size(feature_mat_2,feature_time_ind));

index_here= 1:5:min(size(feature_1_inds,2),size(feature_2_inds,2));

feature_1_inds_rev = 1:numel(index_here);
feature_2_inds_rev = (numel(index_here)+1):2*numel(index_here);

if (plotflag)
figure(999)
plot(SCORE(feature_1_inds(index_here),2),SCORE(feature_1_inds(index_here),1),'b+')
hold on
plot(SCORE(feature_2_inds(index_here),2),SCORE(feature_2_inds(index_here),1),'r+')
hold off


figure(1000)
axis_hist = -200:5:200;
for ind_pick = 1:maxinds
    subplot_tight(4,4,ind_pick)
    [n,x] = hist( SCORE(feature_1_inds,ind_pick),axis_hist);
  %  bar(x,log10(n),'b','EdgeColor','none');
       semilogy(x,log10(n),'b','Linewidth',2);
     hold on
    [n2,x2] = hist( SCORE(feature_2_inds,ind_pick),axis_hist);
   % bar(x2,log10(n2),'r','EdgeColor','none');
        semilogy(x2,log10(n2),'r','Linewidth',2);
    hold off
    if (ind_pick == maxinds)
        legend('Pre','Post')
    end
end
print('-dpng',strcat(plotdir_here,'posture_eigenvalue_hists.png'))
end




gridinds = -200:20:200;
histmat_1  = hist2(SCORE(feature_1_inds(index_here),1),SCORE(feature_1_inds(index_here),2), gridinds, gridinds);
histmat_2  = hist2(SCORE(feature_2_inds(index_here),1),SCORE(feature_2_inds(index_here),2), gridinds, gridinds);

high_counts = intersect(find(histmat_1>10),find(histmat_2>10));
if (plotflag)
figure(1001)
imagesc(histmat_2./histmat_1)
colorbar;
end

ratio_vals = (histmat_2./histmat_1);
[high_zscore,inds] = sort(ratio_vals(high_counts),'DESCEND');

if high_counts
for ind_pick = 1:maxinds;
    if ind_pick <= max(size(inds))
[xhere,yhere] = ind2sub(size(histmat_1),high_counts(inds(ind_pick)));

vals = intersect(intersect(find(SCORE(feature_2_inds(index_here),1)>gridinds(xhere)),find(SCORE(feature_2_inds(index_here),1)<gridinds(xhere+1))),...
    intersect(find(SCORE(feature_2_inds(index_here),2)>gridinds(yhere)),find(SCORE(feature_2_inds(index_here),2)<gridinds(yhere+1))));


vals2 = intersect(intersect(find(SCORE(feature_1_inds(index_here),1)>gridinds(xhere)),find(SCORE(feature_1_inds(index_here),1)<gridinds(xhere+1))),...
    intersect(find(SCORE(feature_1_inds(index_here),2)>gridinds(yhere)),find(SCORE(feature_1_inds(index_here),2)<gridinds(yhere+1))));

if (plotflag)
h = figure (444);
subplot(4,4,ind_pick)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(vals,:),1),3,[])',mod_2.cluster_markersets{clustering_type},mocapstruct_post.markercolor,mocapstruct_post.links,h)
hold on
if (ind_pick==1)
title('no Lesion')
end

hh = figure (445);
subplot(4,4,ind_pick)
plot_eigenpose_subset(reshape(nanmean(feature_mat_2(vals2,:),1),3,[])',mod_2.cluster_markersets{clustering_type},mocapstruct_post.markercolor,mocapstruct_post.links,hh)
hold on
if (ind_pick==1)
title('Lesion')
end
    end
end
end
end

output_struct = [];


%% use other embedding approaches to find the 'most different' poses 
num_dim = 20;
output = cat(2,(ones(1,numel(feature_1_inds(index_here)))),2*(ones(1,numel(feature_1_inds(index_here)))));
Mdl = fitcdiscr(SCORE(cat(2,feature_1_inds(index_here),feature_2_inds(index_here)), find(EXPLAINED>0.1)),output,'Crossval','on','DiscrimType','quadratic');
if (Mdl.KFold>0)
fprintf('Kfold loss %f \n',Mdl.kfoldLoss)
[label,posterior] = kfoldPredict(Mdl);

%Mdl = fitrlinear(SCORE(cat(2,feature_1_inds(index_here),feature_2_inds(index_here)), find(EXPLAINED>0.1)),output,'Crossval','on');
%fprintf('Kfold loss %f \n',Mdl.kfoldLoss)
% [label] = kfoldPredict(Mdl);


if (plotflag)
figure(44)
%[label,posterior] = Mdl.kfoldPredict(Mdl,SCORE(feature_2_inds(index_here),1:10));
plot(posterior)
end
% 
% fprintf('Dimension zscores \n')
% posturechangezscore = abs(Mdl.Trained{1}.Mu(2,:)-Mdl.Trained{1}.Mu(1,:))'./diag(Mdl.Trained{1}.Sigma);
% disp(posturechangezscore)
% [~,ind_zs] = sort(posturechangezscore,'descend');
% 
% posturedelta = Mdl.Trained{1}.Mu(2,:)-Mdl.Trained{1}.Mu(1,:);

ind_zs = 1:numel(find(EXPLAINED>0.1));
posturedelta = zeros(1,numel(find(EXPLAINED>0.1)));
posturechangezscore = zeros(1,numel(find(EXPLAINED>0.1)));

if (plotflag)
h=figure(667)
for ind_pick = 1:4
subplot(1,4,ind_pick)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'w'},num_colors ,1),mocapstruct_post.links,h)
hold on
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])'+posturedelta(ind_zs(ind_pick))*5.*reshape(COEFF(:,ind_zs(ind_pick)),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'r'},num_colors ,1),mocapstruct_post.links,h)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])'-posturedelta(ind_zs(ind_pick))*5.*reshape(COEFF(:,ind_zs(ind_pick)),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'b'},num_colors ,1),mocapstruct_post.links,h)
hold off
ntitle(strcat('explained ',num2str(EXPLAINED(ind_zs(ind_pick)))),'location','north','color','w')
xlabel(strcat('zscore ',num2str(posturechangezscore(ind_zs(ind_pick)))))
view(6,41)
title('red is post')
end
print('-dpng',strcat(plotdir_here,'posture_eigenchange_zscores.png'))
end
output_struct.kfoldloss = Mdl.kfoldLoss;
output_struct.good_post_post =  intersect(find(posterior>0.5),feature_2_inds_rev)-numel(feature_1_inds_rev);
output_struct.good_post_pre =  intersect(find(posterior<0.5),feature_2_inds_rev)-numel(feature_1_inds_rev);
output_struct.good_pre_post =  intersect(find(posterior>0.5),feature_1_inds_rev);
output_struct.good_pre_pre =  intersect(find(posterior<0.5),feature_1_inds_rev);

output_struct.posterior = posterior;
output_struct.feature_1_inds = feature_1_inds;
output_struct.feature_2_inds = feature_2_inds;



[~,framesplot] = sort(posterior(output_struct.good_post_pre),'ascend');

if (plotflag)
h=figure(4444)
pictoral_animation_subset(permute(reshape(feature_mat_2(index_here,:),[],3,numel(mod_2.cluster_markersets{clustering_type})),[1 3 2]),mocapstruct_post,framesplot(1)',50,mod_2.cluster_markersets{clustering_type},h)
end
[~,framesplot] = sort(posterior(output_struct.good_pre_pre),'ascend');

%h=figure(5555)
%pictoral_animation_subset(permute(reshape(feature_mat_1(index_here,:),[],3,numel(mod_2.cluster_markersets{clustering_type})),[1 3 2]),mocapstruct_pre,framesplot(1)',50,mod_2.cluster_markersets{clustering_type},h)

good_post = [];
good_pre = [];
else
    output_struct.kfoldloss = nan;
end
    




elseif strcmp(comptype,'dynamics')

%Mdl2 = fitcsvm(SCORE(cat(2,feature_1_inds(index_here),feature_2_inds(index_here)),1:10),output);

%% do the same analyses for dynamics

opts.fps =300;
opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;
nummarkershere = size(feature_mat_1,2)./3;

agg_features = permute(reshape(feature_mat_1,size(feature_mat_1,1),3,[]),[1 3 2]);
% frames_ind = bsxfun(@plus,1:25:500,floor(2.28*10^5));
% offset = 75;
% 
% h=figure(145);
% pictoral_animation_subset(agg_features,mocapstruct_pre,frames_ind,offset,mod_2.cluster_markersets{clustering_type},h)
% subplot(2,1,1)
% xlim([-100 1500])
% view([-1 16])

totaltime = min(size(feature_mat_1,1),size(feature_mat_2,1));

if (totaltime < 3*opts.fps)
    fprintf('not enough samples \n')
    output_struct.kfoldloss = nan;
  return
end
[dyadic_spectrograms_1,fr,timespect] = get_dyadic_spectrogram(feature_mat_1(1:totaltime,:)',opts);
[dyadic_spectrograms_2,fr,timespect2] = get_dyadic_spectrogram(feature_mat_2(1:totaltime,:)',opts);
feature_time_ind=2;


time_spect_full = cat(2,timespect,timespect2);

feature_1_inds = 1:size(dyadic_spectrograms_1,feature_time_ind);
feature_2_inds = (max(feature_1_inds)+1):(max(feature_1_inds)+size(dyadic_spectrograms_2,feature_time_ind));

if (plotflag)
figure(55)
subplot(2,1,1)
plot(feature_mat_1(:,3))
subplot(2,1,2)
plot(feature_mat_2(:,3))
end

%features from both combined
[COEFF, SCORE_dyn, LATENT, TSQUARED,EXPLAINED] = pca(cat(feature_time_ind,dyadic_spectrograms_1,dyadic_spectrograms_2)');


%% get rid of potential variance explained by shifts in marker position
%[COEFF, SCORE_1, LATENT, TSQUARED,EXPLAINED] = pca(dyadic_spectrograms_1');

%SCORE_2 = bsxfun(@minus,dyadic_spectrograms_1',mean(dyadic_spectrograms_2',1))*COEFF;
%SCORE_dyn = cat(1,SCORE_1,SCORE_2);



if plotflag
maxinds =16;
figure(1000)
axis_hist = -200:5:200;
for ind_pick = 1:maxinds
    subplot_tight(4,4,ind_pick)
    [n,x] = hist( SCORE_dyn(feature_1_inds,ind_pick),axis_hist);
  %  bar(x,log10(n),'b','EdgeColor','none');
       semilogy(x,log10(n),'b','Linewidth',2);
     hold on
    [n2,x2] = hist( SCORE_dyn(feature_2_inds,ind_pick),axis_hist);
   % bar(x2,log10(n2),'r','EdgeColor','none');
        semilogy(x2,log10(n2),'r','Linewidth',2);
    hold off
    if (ind_pick == maxinds)
        legend('Pre','Post')
    end
end
end

if (plotflag)
%disp('variance explained')
%disp(EXPLAINED(1:100))
end
%index_here = 1:min(size(dyadic_spectrograms_1,2),size(dyadic_spectrograms_2,2));

%% use other embedding approaches to find the 'most different' poses 
pcstofit = find(EXPLAINED>0.1);
output = cat(2,(ones(1,numel(feature_1_inds))),2*(ones(1,numel(feature_2_inds))));
Mdl = fitcdiscr(SCORE_dyn(cat(2,feature_1_inds,feature_2_inds),pcstofit),output,'Crossval','on','DiscrimType','quadratic');
fprintf('Kfold loss %f \n',Mdl.kfoldLoss)

fprintf('Dimension zscores \n') 
posturechangezscore = abs(Mdl.Trained{1}.Mu(2,:)-Mdl.Trained{1}.Mu(1,:))'./diag(Mdl.Trained{1}.Sigma);
%disp(posturechangezscore);
[~,ind_zs] = sort(posturechangezscore,'descend');


[label,posterior] = kfoldPredict(Mdl);
if plotflag
figure(49)
%[label,posterior] = Mdl.kfoldPredict(Mdl,SCORE(feature_2_inds(index_here),1:10));
plot(posterior)

figure(45)
for ind_pick = 1:8
subplot_tight(2,4,ind_pick)
%imagesc(reshape(COEFF(:,ind_zs(ind_pick)),[],36))
h=pcolor(1:nummarkershere*3, fr, reshape(COEFF(:,ind_zs(ind_pick)),[],nummarkershere*3)); 
set(h, 'EdgeColor', 'none');
end



figure(460)
for ind_pick = 1:8
subplot(2,4,ind_pick)
reshaped_vals = reshape(COEFF(:,ind_zs(ind_pick)),[],nummarkershere*3);
reshaped_vals(reshaped_vals<0) = nan;
plot(fr,nansum(reshaped_vals,2))

end


figure(470)
for ind_pick = 1:8
subplot(2,4,ind_pick)
reshaped_vals_here = sum(reshape(COEFF(:,ind_zs(ind_pick)),[],nummarkershere*3),1);
bar(reshaped_vals_here)
set(gca,'XTick',1:3:nummarkershere*3,'XTickLabel',mocapstruct_pre.markernames(mod_1.cluster_markersets{clustering_type}) )

end

figure(480)
for ind_pick = 1:8
subplot(2,4,ind_pick)
reshaped_vals_here = sum(reshape(COEFF(:,ind_zs(ind_pick)),[],nummarkershere*3),1);
xyz_inds = cat(2,1:3:numel(reshaped_vals_here),2:3:numel(reshaped_vals_here),3:3:numel(reshaped_vals_here));
bar(reshaped_vals_here(xyz_inds))
set(gca,'XTick',1:1:nummarkershere,'XTickLabel',mocapstruct_pre.markernames(mod_1.cluster_markersets{clustering_type}) )

end
end

%get top 100 values of SCORE and look at videos around these timepoints
%[vals,inds] = sort(SCORE(:,ind_zs(1)),'DESCEND');

[vals,inds] = sort(posterior(:,2),'DESCEND');

maxind = floor(size(posterior,1)./2);
inds_pick_post = inds(1:maxind ,1);
inds_pick_pre = inds(size(inds)-maxind :end,1);

output_struct = [];
output_struct.good_post_post =  frame_ex_2(time_spect_full(intersect(find(output == 2),inds_pick_post))*opts.fps);
output_struct.good_post_pre =  frame_ex_1(time_spect_full(intersect(find(output == 1),inds_pick_post))*opts.fps);

output_struct.good_pre_post =  frame_ex_2(time_spect_full(intersect(find(output == 2),inds_pick_pre))*opts.fps);
output_struct.good_pre_pre =  frame_ex_1(time_spect_full(intersect(find(output == 1),inds_pick_pre))*opts.fps);
output_struct.kfoldloss = Mdl.kfoldLoss;

% good_post =output_struct.good_post_post;
% good_pre =output_struct.good_pre_pre;

end

end


%% train a linear discriminant and look at the prediction for differences and 
