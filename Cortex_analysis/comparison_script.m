
%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
createmocapfilestruct('JDM25',mocapmasterdirectory) %this step can take an hour, potentially longer on the server

createmocapfilestruct('JDM32',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct_JDM25 = loadmocapfilestruct('JDM25',mocapmasterdirectory);
mocapfilestruct_JDM32 = loadmocapfilestruct('JDM32',mocapmasterdirectory);

mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);

[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_prelesion',mocapmasterdirectory);
[mocapstruct_pre] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,0,0);
mocapstruct_pre = mocapstruct_pre.mocap_struct;

[descriptor_struct_2,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_dlslesion_late',mocapmasterdirectory);
[mocapstruct_post] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes,0,1);


[descriptor_struct_3,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25bipost',mocapmasterdirectory);
[mocapstruct_post_bi] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes,0,1);


plotfractionmissing(mocapstruct_pre)

compare_plot_marker_characteristics_timerange(mocapstruct_pre,mocapstruct_pre.move_frames,mocapstruct_post,mocapstruct_post.move_frames)

compare_comp_fraction_moving(mocapstruct_pre,mocapstruct_post_bi)
%compare_plot_marker_characteristics(mocapstruct)

%compare the eigen clusters based on specific timepoins
cluster_here = [8];
[modular_cluster_properties] = get_modularclustproperties(mocapstruct_pre);
[modular_cluster_properties2] = get_modularclustproperties(mocapstruct_post);
[modular_cluster_properties3] = get_modularclustproperties(mocapstruct_post_bi);

predict_markers_position(mocapstruct_pre,modular_cluster_properties)
predict_markers_position(mocapstruct_post,modular_cluster_properties2)
predict_markers_position(mocapstruct_post_bi,modular_cluster_properties3)


[output_struct] = compare_features(mocapstruct_pre,modular_cluster_properties,mocapstruct_post,modular_cluster_properties2,'dynamics');

pre_inds = 1:numel(frame_ex_1);
post_inds = (numel(frame_ex_1)+1):(numel(frame_ex_1)+numel(frame_ex_2));


good_post_sort = sort(good_post,'ASCEND');
good_post_sort = unique(bsxfun(@plus,good_post_sort,(-150:150)'));


good_pre_sort = sort(good_pre,'ASCEND');
good_pre_sort = unique(bsxfun(@plus,good_pre_sort,(-150:150)'));

% 
% 
% good_post_post = intersect(good_post_sort,post_inds);
% good_post_pre = intersect(good_post_sort,pre_inds);
% 
% 
% good_pre_post = intersect(good_pre_sort,post_inds);
% good_pre_pre = intersect(good_pre_sort,pre_inds);

 dos('cd C:\Users\Jesse Marshall\Documents\GitHub\Movement_analysis\FATGUI-master\')
MCC;

%% the clipped index is the index of clipped frames in the full trace
[modular_cluster_properties] = get_clustering_features(mocapstruct_pre,modular_cluster_properties,cluster_here)  ;
[modular_cluster_properties2] = get_clustering_features(mocapstruct_post,modular_cluster_properties2,cluster_here)  ;

[licktime,levertime] = gettasktime(mocapstruct);
leverind = find(licktime==1);

animate_markers_aligned_fullmovie(mocapstruct_post,modular_cluster_properties2.clustering_inds_agg{2}(1:20:end))

%% get the spectral PCS of each and embed? 

    %% first compute the eigenposes to look for postural differences
      %% look at the eigenposes of the data
        fprintf('computing eigenposes \n')
  time_subset = mocapstruct_pre.move_frames;
  time_subset2 = mocapstruct_post.move_frames;
  
  %high rear
  time_subset = find(( mocapstruct_pre.markers_preproc.SpineF(:,3)-mocapstruct_pre.markers_preproc.SpineL(:,3))>50);
  time_subset2 = find(( mocapstruct_post.markers_preproc.SpineF(:,3)-mocapstruct_post.markers_preproc.SpineL(:,3))>50);
  
  %low rear -- shortens stance more and more
  time_subset = find(( mocapstruct_pre.markers_preproc.HeadB(:,3)-mocapstruct_pre.markers_preproc.SpineF(:,3))<-50);
  time_subset2 = find(( mocapstruct_post.markers_preproc.HeadB(:,3)-mocapstruct_post.markers_preproc.SpineF(:,3))<-50);
  
  %l/r groom 
  time_subset = find(( mocapstruct_pre.markers_aligned_preproc.SpineF(:,2)-mocapstruct_pre.markers_aligned_preproc.SpineL(:,2))>20);
  time_subset2 = find(( mocapstruct_post.markers_aligned_preproc.SpineF(:,2)-mocapstruct_post.markers_aligned_preproc.SpineL(:,2))>20);
  
  %look at distributions of distances
  
   % time_subset = find(( markers_clipped.SpineL(:,2)>20));
  
 eigenpose = compare_eigenposes(mocapstruct_pre,modular_cluster_properties,cluster_here,time_subset,...
     mocapstruct_post,modular_cluster_properties2,time_subset2);
 
 
 % time_subset = mocapstruct_pre.move_frames;
 % time_subset2 = mocapstruct_post.move_frames;
 
  %% break this up into its own script
 figure(33)
%[n,x] = hist(( mocapstruct_pre.markers_aligned_preproc.HeadF(:,2)-mocapstruct_pre.markers_aligned_preproc.SpineF(:,2)),-100:1:100);
%[n2,x2] = hist(( mocapstruct_post.markers_aligned_preproc.HeadF(:,2)-mocapstruct_post.markers_aligned_preproc.SpineF(:,2)),-100:1:100);
%[n,x] = hist(( mocapstruct_pre.markers_aligned_preproc.SpineF(time_subset,3)-mocapstruct_pre.markers_aligned_preproc.SpineM(time_subset,3)),-100:1:100);
%[n2,x2] = hist(( mocapstruct_post.markers_aligned_preproc.SpineF(time_subset2,3)-mocapstruct_post.markers_aligned_preproc.SpineM(time_subset2,3)),-100:1:100);
%[n,x] = hist(( mocapstruct_pre.markers_aligned_preproc.SpineM(time_subset,2)-mocapstruct_pre.markers_aligned_preproc.SpineL(time_subset,2)),-100:1:100);
%[n2,x2] = hist(( mocapstruct_post.markers_aligned_preproc.SpineM(time_subset2,2)-mocapstruct_post.markers_aligned_preproc.SpineL(time_subset2,2)),-100:1:100);

[n,x] = hist(( mocapstruct_pre.markers_aligned_preproc.KneeL(time_subset,1)-mocapstruct_pre.markers_aligned_preproc.KneeR(time_subset,1)),-100:1:100);
[n2,x2] = hist(( mocapstruct_post.markers_aligned_preproc.KneeL(time_subset2,1)-mocapstruct_post.markers_aligned_preproc.KneeR(time_subset2,1)),-100:1:100);


[n,x] = hist(( mocapstruct_pre.markers_aligned_preproc.ShoulderL(time_subset,3)-mocapstruct_pre.markers_aligned_preproc.ShoulderR(time_subset,3)),-100:1:100);
[n2,x2] = hist(( mocapstruct_post.markers_aligned_preproc.ShoulderL(time_subset2,3)-mocapstruct_post.markers_aligned_preproc.ShoulderR(time_subset2,3)),-100:1:100);
bar(x,n./sum(n),'EdgeColor','none')
hold on
b1 = bar(x2,n2./sum(n2),'r','EdgeColor','none')
b1.FaceAlpha = 0.5;

hold off


compare_plot_marker_characteristics_timerange(mocapstruct_pre,time_subset,mocapstruct_post,time_subset2)


% 
% 
% annotation_struct = struct();
% make_arena_line_plots(mocapstruct)
%make_circadian_plots(mocapstruct,annotation_struct)