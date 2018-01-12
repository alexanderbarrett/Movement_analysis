
%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
createmocapfilestruct('JDM25',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
%createmocapfilestruct('JDM21',mocapmasterdirectory) 
mocapfilestruct_JDM21 = loadmocapfilestruct('JDM21',mocapmasterdirectory);


createmocapfilestruct('JDM32',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct_JDM25 = loadmocapfilestruct('JDM25',mocapmasterdirectory);
mocapfilestruct_JDM32 = loadmocapfilestruct('JDM32',mocapmasterdirectory);
mocapfilestruct_Vicon8 = loadmocapfilestruct('Vicon8',mocapmasterdirectory);


%preprocessallmocap('JDM32',mocapmasterdirectory)
%preprocessallmocap('JDM25',mocapmasterdirectory)


mocapstructcell = cell(1,10);
%% load desired days for comparison

ov_meta = 0;

descriptor_struct = get_mocap_files_table(6,'Vicon8');
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_v8pl2] = preprocess_mocap_data_2(mocapfilearray(5),mocapfilestruct,descriptor_struct,mocapfiletimes,1,0,[],mocapvideodirectory,0);
 base = 300000;
M= animate_markers_aligned_fullmovie_syncedvideo(mocapstruct_v8pl2, base+(1:10:10000),...
    mocapstruct_v8pl2.cameradirectory{2},mocapstruct_v8pl2.matched_frames_aligned{2}( base+(1:10:10000)));

%vicon8 
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM21','JDM21',mocapmasterdirectory);
[mocapstruct_21] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,0,1);

%vicon8 
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_prelesion',mocapmasterdirectory);
[mocapstructcell{9}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,0,0);

[descriptor_struct_2,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_dlslesion_late',mocapmasterdirectory);
[mocapstructcell{10}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes,0,ov_meta);

[descriptor_struct_2,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_dlslesion_late2',mocapmasterdirectory);
[mocapstructcell{16}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes,0,ov_meta);
%caff
[descriptor_struct_2,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);
[mocapstructcell{11}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes,0,ov_meta);
%amph
[descriptor_struct_2,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_amph',mocapmasterdirectory);
[mocapstructcell{12}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes,0,ov_meta);

%second prelesion
[descriptor_struct_2,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_prelesion2',mocapmasterdirectory);
[mocapstructcell{13}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes,0,ov_meta);




%% jdm 25
[descriptor_struct_4,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25pre',mocapmasterdirectory);
[mocapstructcell{4}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_4,mocapfiletimes,0,0);
%mocapstruct_prejdm = mocapstruct_prejdm.mocap_struct;

%unilateral
[descriptor_struct_5,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25post',mocapmasterdirectory);
[mocapstructcell{5}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_5,mocapfiletimes,0,ov_meta);
%mocapstruct_post_uni = mocapstruct_post_uni.mocap_struct;
%bilateral
[descriptor_struct_3,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25bipost',mocapmasterdirectory);
[mocapstructcell{6}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_3,mocapfiletimes,0,ov_meta);
%mocapstruct_post_bi = mocapstruct_post_bi.mocap_struct;

[descriptor_struct_3,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25bipost2',mocapmasterdirectory);
[mocapstructcell{15}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_3,mocapfiletimes,0,ov_meta);

%pre2
[descriptor_struct_4,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25pre2',mocapmasterdirectory);
[mocapstructcell{14}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_4,mocapfiletimes,0,1);

[descriptor_struct_3,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25bipost_task',mocapmasterdirectory);
[mocapstructcell{21}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_3,mocapfiletimes,0,ov_meta);

[descriptor_struct_4,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25pre_task',mocapmasterdirectory);
[mocapstructcell{22}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_4,mocapfiletimes,0,ov_meta);

%all

[descriptor_struct_3,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25bipost_all',mocapmasterdirectory);
[mocapstructcell{25}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_3,mocapfiletimes,0,ov_meta);

[descriptor_struct_4,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25pre_all',mocapmasterdirectory);
[mocapstructcell{26}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_4,mocapfiletimes,0,ov_meta);

%caff pre
[descriptor_struct_t,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25caff',mocapmasterdirectory);
[mocapstructcell{7}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_t,mocapfiletimes,0,0);
%mocapstruct_caff_jdm25=mocapstruct_caff_jdm25.mocap_struct;
%amph pre
[descriptor_struct_t,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25amph',mocapmasterdirectory);
[mocapstructcell{8}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_t,mocapfiletimes,0,ov_meta);
%mocapstruct_amph_jdm25=mocapstruct_amph_jdm25.mocap_struct;

%pre long
%[descriptor_struct_t,mocapfilearray,mocapfilestruct_JDM25,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM25','JDM25pre_long',mocapmasterdirectory);
%[mocapstruct_pre_long_jdm25] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM25,descriptor_struct_t,mocapfiletimes,0,ov_meta);



%% jdm 32
[descriptor_struct_6,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32postbi_2',mocapmasterdirectory);
[mocapstructcell{1}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_6,mocapfiletimes,0,ov_meta);

[descriptor_struct_7,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32pre',mocapmasterdirectory);
[mocapstructcell{2}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_7,mocapfiletimes,0,ov_meta);

%uni
[descriptor_struct_8,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32postuni',mocapmasterdirectory);
[mocapstructcell{3}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_8,mocapfiletimes,0,ov_meta);
%caff
%amph
[descriptor_struct_7,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32postbi_3',mocapmasterdirectory);
[mocapstructcell{17}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_7,mocapfiletimes,0,ov_meta);

[descriptor_struct_7,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32pre2',mocapmasterdirectory);
[mocapstructcell{18}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_7,mocapfiletimes,0,ov_meta);

%task
[descriptor_struct_7,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32postbi_task',mocapmasterdirectory);
[mocapstructcell{19}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_7,mocapfiletimes,0,ov_meta);

[descriptor_struct_7,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32pre_task',mocapmasterdirectory);
[mocapstructcell{20}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_7,mocapfiletimes,0,ov_meta);

%all
[descriptor_struct_7,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32postbi_all',mocapmasterdirectory);
[mocapstructcell{23}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_7,mocapfiletimes,0,ov_meta);

[descriptor_struct_7,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32pre_all',mocapmasterdirectory);
[mocapstructcell{24}] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_7,mocapfiletimes,0,0);



mocapstruct_types = {'jdm32_bi','jdm32_pre','jdm32_uni',...
    'jdm25_pre','jdm25_uni','jdm25_postbi','jdm25_caff','jdm25_amph','vicon8_pre','vicon8_uni','vicon8_caff',...
    'vicon8_amph','vicon8_pre2','jdm25_pre2','jdm25_bipost2','vicon8_uni2','jdm32_post2','jdm32pre2',...
    'jdm32_bitask','jdm32_pretask','jdm25_bitask','jdm25_pretask','jdm32_biall','jdm32_preall','jdm25_biall','jdm25_preall'};


comp_list_32 = [1,2,3,17,18];%[1:3];
comp_list_25 = [4,5,6,14,15];%[4:8,14];
comp_list_8 = [9,10,13,16];%9:13;
comp_list_task_25 = [21,22];
comp_list_task_32 = [19,20];
comp_list_all_25 = [25,26];
comp_list_all_32 = [23,24];
comp_list_25_mix = [4,14,22,6,15,21];%[4:8,14];
comp_list_32_mix = [2,18,20,1,17,19];%[4:8,14];

%% get the desired recordings

%% compare gross characteristics between conditions
plotfractionmissing(mocapstruct_pre)

compare_plot_marker_characteristics_timerange(mocapstruct_pre,mocapstruct_pre.move_frames,mocapstruct_post,mocapstruct_post.move_frames)


compare_plot_marker_characteristics_timerange(mocapstruct_prejdm,mocapstruct_prejdm.move_frames,mocapstruct_post_bi,mocapstruct_post_bi.move_frames)

compare_comp_fraction_moving(mocapstruct_prejdm,mocapstruct_post_uni32)
compare_comp_fraction_moving(mocapstruct_prejdm,mocapstruct_post_uni)

colorarray = {'r','b','g','k','c','y'};

list_compare = comp_list_32_mix
legend_tags = mocapstruct_types(list_compare);
for kk =list_compare
    kk_ind = find(list_compare == kk);
multi_comp_fraction_moving(mocapstructcell{kk},colorarray{kk_ind},kk_ind)
%timestruct{kk} = mocapstructcell{kk}.move_frames;
%multi_plot_marker_characteristics_timerange(mocapstructcell{kk},timestruct{kk},colorarray{kk_ind},kk_ind)

end

figure(105)
subplot(2,1,1)
legend(legend_tags)
subplot(2,1,2)
legend(legend_tags)
% save here
figure(306)
legend(legend_tags)
figure(307)
subplot(1,2,1)
legend(legend_tags)
subplot(1,2,2)
legend(legend_tags)

 %% for the PSD comparison/multi plot characteristics
 figure(200)
 legend(legend_tags)
 
  figure(300)
 legend(legend_tags)
 
%% Look at the outliers 
 spinevel = get_spine_vel(mocapstructcell{4});
 velinds = intersect(find(spinevel>5),find(spinevel<7));
 animate_markers_aligned_fullmovie(mocapstructcell{4},velinds(1:10:end));
 
  spinevel = get_spine_vel(mocapstructcell{2});
 velinds = intersect(find(spinevel>15),find(spinevel<25));
 animate_markers_aligned_fullmovie(mocapstructcell{2},velinds(1:10:end));
 
% 
 animate_markers_aligned_fullmovie(mocapstruct_21,8*10^5:10:10*10^5);

mocapstruct_21

cluster_here = [2];

mocapstruct_types = {'jdm32_bi','jdm32_pre','jdm32_uni',...
    'jdm25_pre','jdm25_uni','jdm25_postbi','jdm25_caff','jdm25_amph','vicon8_pre','vicon8_uni','vicon8_caff','vicon8_amph'};


comp_list_32 = [1:3];
comp_list_25 = 4:8;

%% get the spectral PCS of each and embed? 

    %% first compute the eigenposes to look for postural differences
      %% look at the eigenposes of the data
        fprintf('computing eigenposes \n')

    %compare the frequencies
  compare_pose_frequencies(mocapstructcell{4}.annotated_pose_struct,mocapstructcell{5}.annotated_pose_struct, mocapstructcell{6}.annotated_pose_struct,...
      mocapstructcell{14}.annotated_pose_struct,mocapstructcell{15}.annotated_pose_struct,...
      {'pre','uni','bi post','pre2','bipost2'})
  
%    compare_pose_frequencies(mocapstructcell{9}.annotated_pose_struct,mocapstructcell{10}.annotated_pose_struct, mocapstructcell{11}.annotated_pose_struct,...
%       mocapstructcell{12}.annotated_pose_struct,...
%       {'pre','uni','caff','amph'})
%   
   compare_pose_frequencies(mocapstructcell{9}.annotated_pose_struct,mocapstructcell{10}.annotated_pose_struct, mocapstructcell{13}.annotated_pose_struct,...
      mocapstructcell{16}.annotated_pose_struct,...
      {'pre','uni','pre2','uni2'})
  
   compare_pose_frequencies(mocapstructcell{1}.annotated_pose_struct,mocapstructcell{2}.annotated_pose_struct, mocapstructcell{3}.annotated_pose_struct,...
      mocapstructcell{17}.annotated_pose_struct,mocapstructcell{18}.annotated_pose_struct,...
      {'bi','pre','uni','bi2','pre2'})



%% compare eigenposes

%   eigenpose = compare_eigenposes(mocapstruct_prejdm,modular_cluster_properties4,cluster_here,    annotated_pose_struct_pre.(comparefield),...
%      mocapstruct_post_uni,modular_cluster_properties5,annotated_pose_struct_uni.(comparefield) );
%  comparefield = 'high_rear';
 eigenpose = compare_eigenposes(mocapstructcell{4},mocapstructcell{4}.modular_cluster_properties,2,mocapstructcell{4}.annotated_pose_struct.(comparefield),...
    mocapstructcell{6},mocapstructcell{6}.modular_cluster_properties,mocapstructcell{6}.annotated_pose_struct.(comparefield) );
comparefield = 'high_rear';
 
[output_struct] = compare_features(mocapstructcell{9},mocapstructcell{9}.modular_cluster_properties,...
    mocapstructcell{13},mocapstructcell{13}.modular_cluster_properties,'dynamics',....
   mocapstructcell{9}.annotated_pose_struct.(comparefield),mocapstructcell{13}.annotated_pose_struct.(comparefield),2,1 );

 
 
 

 % time_subset = mocapstruct_pre.move_frames;
 % time_subset2 = mocapstruct_post.move_frames;
 %% loop over and plot distb
 
colorarray = {'r','b','g','k','c','y'};
legend_tags = {'pre','post uni','post bi','caff','amph','pre long'};
cellstruct = cell(1,5);
cellstruct{1} = mocapstruct_post_bi;
cellstruct{2} = mocapstruct_prejdm;
annotstruct{1} = annotated_pose_struct_post.(comparefield);
annotstruct{2} = annotated_pose_struct_pre.(comparefield);

for kk =1:2
 plot_intermarker_distribution(cellstruct{kk},annotstruct{kk},colorarray{kk})
end
 

%% something else

 params.fps = 300;
 pose_fragments_pre = hipass_clip_fragments(mocapstruct_prejdm.markers_aligned_preproc,annotated_pose_struct_pre.(comparefield),params);
  pose_fragments_post = hipass_clip_fragments(mocapstruct_post_bi.markers_aligned_preproc,annotated_pose_struct_post.(comparefield),params);

  figure(33)
  subplot(2,1,1)
plot( pose_fragments_pre.HeadF(:,3))
subplot(2,1,2)
plot( pose_fragments_post.HeadF(:,3))


%% look at predictions to see how modular behavior is
predict_markers_position_split(mocapcellstruct{1},modular_cluster_properties4,'pose')
predict_markers_position(mocapstruct_post,modular_cluster_properties2)
predict_markers_position_split(mocapstruct_post_bi,modular_cluster_properties3,'pose')
predict_markers_position_split(mocapstructcell{22},mocapstructcell{22}.modular_cluster_properties,'dynamics')



%% see if a classifier can differentiate the two and if so, what is it using and which markers can differentiate (POSE)

compare_multi_features(mocapstructcell([9,10,13,16]),2,mocapstruct_types([9,10,13,16]),'pose') %,'pre long','pre 32','caff','amph'
compare_multi_features(mocapstructcell([4,14,5,6,15]),2,mocapstruct_types([4,14,5,6,15]),'pose') %,'pre long','pre 32','caff','amph'
compare_multi_features(mocapstructcell([comp_list_32([2,5,3,1,4])]),2,mocapstruct_types([comp_list_32([2,5,3,1,4])]),'pose') %,'pre long','pre 32','caff','amph'
compare_multi_features(mocapstructcell([20 22]),2,mocapstruct_types([20 22]),'pose') %,'pre long','pre 32','caff','amph'

compare_multi_features(mocapstructcell([23 24]),2,mocapstruct_types([23 24]),'pose') %,'pre long','pre 32','caff','amph'
compare_multi_features(mocapstructcell([comp_list_32_mix]),2,mocapstruct_types([comp_list_32_mix]),'pose') %,'pre long','pre 32','caff','amph'
compare_multi_features(mocapstructcell([comp_list_25_mix]),2,mocapstruct_types([comp_list_25_mix]),'pose') %,'pre long','pre 32','caff','amph'


comp_list_task_25

%% attempt to classify dynamics

comparefield = 'moving';
mocapcellstruct = cell(1,2);
mocapcellstruct{1} = mocapstruct_prejdm32;
mocapcellstruct{2} = mocapstruct_post_bi32;

compare_multi_features(mocapcellstruct,2,{'pre','post'},'dynamics')
animate_markers_aligned_fullmovie(mocapstruct_post_bi32,mocapstruct_post_bi32.annotated_pose_struct.low_rear(1:10:end))
animate_markers_aligned_fullmovie(mocapstruct_prejdm32,intersect(mocapstruct_prejdm32.modular_cluster_properties.clustering_inds_agg{2},....
    mocapstruct_prejdm32.annotated_pose_struct.high_rear(1:10:end)))

%% other single pairwise comparisons

[output_struct] = compare_features(mocapstruct_prejdm,mocapstruct_prejdm.modular_cluster_properties,...
    mocapstruct_post_bi,mocapstruct_post_bi.modular_cluster_properties,'dynamics',....
   mocapstruct_prejdm.annotated_pose_struct.(comparefield),mocapstruct_post_bi.annotated_pose_struct.(comparefield),7,1 );

[output_struct] = compare_features(mocapstruct_prejdm32,mocapstruct_prejdm32.modular_cluster_properties,...
    mocapstruct_post_bi32,mocapstruct_post_bi32.modular_cluster_properties,'dynamics',....
  mocapstruct_prejdm32.annotated_pose_struct.(comparefield),mocapstruct_post_bi32.annotated_pose_struct.(comparefield),7,1 );


[output_struct] = compare_features(mocapstruct_prejdm32,mocapstruct_prejdm32.modular_cluster_properties,...
    mocapstruct_prejdm,mocapstruct_prejdm.modular_cluster_properties,'dynamics',....
  mocapstruct_prejdm32.annotated_pose_struct.(comparefield),mocapstruct_prejdm.annotated_pose_struct.(comparefield),7,1 );



[output_struct] = compare_features(mocapstruct_prejdm32,modular_cluster_properties6,mocapstruct_post_bi32,modular_cluster_properties7,'dynamics',...
    annotated_pose_struct_pre32.(comparefield),annotated_pose_struct_post32.(comparefield));



%% visualize frames that I claim are different

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


% 
% 
% annotation_struct = struct();
% make_arena_line_plots(mocapstruct)
%make_circadian_plots(mocapstruct,annotation_struct)