function [good_post,good_pre,frame_ex_1,frame_ex_2] = compare_features(mocapstruct_pre,mod_1,mocapstruct_post,mod_2,comptype)

    struct_fn = fieldnames(mocapstruct_pre.markers_aligned_preproc);
    clustering_type =2;
    feature_mat_1 = [];
        feature_mat_2 = [];

        plotdir_here = strcat(mocapstruct_pre.plotdirectory,filesep,'comparison_subcluster',num2str(clustering_type),filesep);
        mkdir(plotdir_here);
        
    downsample = 3;

    for mm = mod_1.cluster_markersets{clustering_type}    
        
 frame_ex_1 = intersect(mod_1.clustering_inds_agg{clustering_type}(1:downsample:end),mocapstruct_pre.move_frames);
    frame_ex_2 = intersect(mod_2.clustering_inds_agg{clustering_type}(1:downsample:end),mocapstruct_post.move_frames);
     frames_total = cat(2,frame_ex_1,frame_ex_2);
     
     feature_mat_1 = cat(2,feature_mat_1,mocapstruct_pre.markers_aligned_preproc.(struct_fn{mm})(frame_ex_1,:));
          feature_mat_2 = cat(2,feature_mat_2,mocapstruct_post.markers_aligned_preproc.(struct_fn{mm})(frame_ex_2,:));
    end
             % feature_mat_2 = feature_mat_1(250000:end,:);

    % feature_mat_1 = feature_mat_1(1:250000,:);
if (strcmp(comptype,'pose'))


%% first compute PCS of the joint and look at sample 2 and 3d point clouds
feature_time_ind = 1;
[COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(cat(feature_time_ind,feature_mat_1,feature_mat_2));

maxinds = min(size(COEFF,1),15);

%% plot the PCs
num_colors = numel(mocapstruct_post.markercolor);

h=figure(666)
for ind_pick = 1:maxinds
subplot(4,4,ind_pick)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'w'},num_colors ,1),mocapstruct_post.links,h)
hold on
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])'+50.*reshape(COEFF(:,ind_pick),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'r'},num_colors ,1),mocapstruct_post.links,h)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(:,:),1),3,[])'-50.*reshape(COEFF(:,ind_pick),3,[])',mod_2.cluster_markersets{clustering_type},repmat({'b'},num_colors ,1),mocapstruct_post.links,h)

hold off
end




feature_1_inds = 1:size(feature_mat_1,feature_time_ind);
feature_2_inds = (max(feature_1_inds)+1):(max(feature_1_inds)+size(feature_mat_2,feature_time_ind));

index_here= 1:5:min(size(feature_1_inds,2),size(feature_2_inds,2));

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




gridinds = -200:20:200;
histmat_1  = hist2(SCORE(feature_1_inds(index_here),1),SCORE(feature_1_inds(index_here),2), gridinds, gridinds);
histmat_2  = hist2(SCORE(feature_2_inds(index_here),1),SCORE(feature_2_inds(index_here),2), gridinds, gridinds);

high_counts = intersect(find(histmat_1>10),find(histmat_2>10));
figure(1001)
imagesc(histmat_2./histmat_1)
colorbar;

ratio_vals = (histmat_2./histmat_1);
[high_zscore,inds] = sort(ratio_vals(high_counts),'DESCEND');

for ind_pick = 1:maxinds;
[xhere,yhere] = ind2sub(size(histmat_1),high_counts(inds(ind_pick)));

vals = intersect(intersect(find(SCORE(feature_2_inds(index_here),1)>gridinds(xhere)),find(SCORE(feature_2_inds(index_here),1)<gridinds(xhere+1))),...
    intersect(find(SCORE(feature_2_inds(index_here),2)>gridinds(yhere)),find(SCORE(feature_2_inds(index_here),2)<gridinds(yhere+1))));


vals2 = intersect(intersect(find(SCORE(feature_1_inds(index_here),1)>gridinds(xhere)),find(SCORE(feature_1_inds(index_here),1)<gridinds(xhere+1))),...
    intersect(find(SCORE(feature_1_inds(index_here),2)>gridinds(yhere)),find(SCORE(feature_1_inds(index_here),2)<gridinds(yhere+1))));

h = figure (444);
subplot(4,4,ind_pick)
plot_eigenpose_subset(reshape(nanmean(feature_mat_1(vals,:),1),3,[])',mod_2.cluster_markersets{clustering_type},mocapstruct_post.markercolor,mocapstruct_post.links,h)
hold on

hh = figure (445);
subplot(4,4,ind_pick)
plot_eigenpose_subset(reshape(nanmean(feature_mat_2(vals2,:),1),3,[])',mod_2.cluster_markersets{clustering_type},mocapstruct_post.markercolor,mocapstruct_post.links,hh)
hold on
end



%% use other embedding approaches to find the 'most different' poses 
output = cat(2,(ones(1,numel(feature_1_inds(index_here)))),2*(ones(1,numel(feature_1_inds(index_here)))));
Mdl = fitcdiscr(SCORE(cat(2,feature_1_inds(index_here),feature_2_inds(index_here)),1:15),output,'Crossval','on','DiscrimType','linear');
fprintf('Kfold loss %f \n',Mdl.kfoldLoss)

figure(44)
[label,posterior] = kfoldPredict(Mdl);
%[label,posterior] = Mdl.kfoldPredict(Mdl,SCORE(feature_2_inds(index_here),1:10));
plot(posterior)


fprintf('Dimension zscores \n') 
posturechangezscore = abs(Mdl.Trained{1}.Mu(2,:)-Mdl.Trained{1}.Mu(1,:))'./diag(Mdl.Trained{1}.Sigma);
disp(posturechangezscore)
[~,ind_zs] = sort(posturechangezscore,'descend');

posturedelta = Mdl.Trained{1}.Mu(2,:)-Mdl.Trained{1}.Mu(1,:);


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
end
print('-dpng',strcat(plotdir_here,'posture_eigenchange_zscores.png'))

good_post = [];
good_pre = [];

elseif strcmp(comptype,'dynamics')

%Mdl2 = fitcsvm(SCORE(cat(2,feature_1_inds(index_here),feature_2_inds(index_here)),1:10),output);

%% do the same analyses for dynamics

opts.fps =300./downsample;
opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;
nummarkershere = size(feature_mat_1,2)./3;


[dyadic_spectrograms_1,fr,timespect] = get_dyadic_spectrogram(feature_mat_1(:,:)',opts);
[dyadic_spectrograms_2,fr,timespect2] = get_dyadic_spectrogram(feature_mat_2(:,:)',opts);
feature_time_ind=2;

time_spect_full = cat(2,timespect,timespect2+max(timespect)+timespect(2)-timespect(1));

feature_1_inds = 1:size(dyadic_spectrograms_1,feature_time_ind);
feature_2_inds = (max(feature_1_inds)+1):(max(feature_1_inds)+size(dyadic_spectrograms_2,feature_time_ind));


[COEFF, SCORE_dyn, LATENT, TSQUARED,EXPLAINED] = pca(cat(feature_time_ind,dyadic_spectrograms_1,dyadic_spectrograms_2)');

index_here = 1:min(size(dyadic_spectrograms_1,2),size(dyadic_spectrograms_2,2));

%% use other embedding approaches to find the 'most different' poses 
pcstofit = 1:15;
output = cat(2,(ones(1,numel(feature_1_inds(index_here)))),2*(ones(1,numel(feature_1_inds(index_here)))));
Mdl = fitcdiscr(SCORE_dyn(cat(2,feature_1_inds(index_here),feature_2_inds(index_here)),pcstofit),output,'Crossval','on','DiscrimType','linear');
fprintf('Kfold loss %f \n',Mdl.kfoldLoss)

fprintf('Dimension zscores \n') 
posturechangezscore = abs(Mdl.Trained{1}.Mu(2,:)-Mdl.Trained{1}.Mu(1,:))'./diag(Mdl.Trained{1}.Sigma);
disp(posturechangezscore)
[~,ind_zs] = sort(posturechangezscore,'descend');


figure(49)
[label,posterior] = kfoldPredict(Mdl);
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


%get top 100 values of SCORE and look at videos around these timepoints
%[vals,inds] = sort(SCORE(:,ind_zs(1)),'DESCEND');


[vals,inds] = sort(posterior(:,2),'DESCEND');
good_post =frames_total(time_spect_full(inds(1:1000,1))*opts.fps);
good_pre =frames_total(time_spect_full(inds(size(inds)-1000:end,1))*opts.fps);

end

end


%% train a linear discriminant and look at the prediction for differences and 
