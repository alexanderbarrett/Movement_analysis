%% annotation test

%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
%mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
loadmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server

mocapfilestruct = loadmocapfilestruct('JDM25',mocapmasterdirectory);

%alternate: Vicon8_caff, Vicon8_dlslesion
[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);
[mocapstruct_caff] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);
ML_features = get_supervised_features(mocapstruct_caff,mocapstruct_caff.modular_cluster_properties.clustering_inds_agg{2},2);
plot_feature_spaces(ML_features)

[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('Vicon8','Vicon8_amph',mocapmasterdirectory);
[mocapstruct_amph] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);
ML_features_amph = get_supervised_features(mocapstruct_amph,mocapstruct_amph.modular_cluster_properties.clustering_inds_agg{2},2);


descriptor_struct = get_mocap_files_table(9,'Vicon8');
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_lesion] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,0);
 %[mocapstruct_lesion] = preprocess_mocap_data2( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,1,[],mocapvideodirectory,0);
ML_features_lesion = get_supervised_features(mocapstruct_lesion,mocapstruct_lesion.modular_cluster_properties.clustering_inds_agg{2},2);



[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('Vicon8','Vicon8_prelesion',mocapmasterdirectory);
[mocapstruct_prelesion] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);
ML_features_prel = get_supervised_features(mocapstruct_prelesion,mocapstruct_prelesion.modular_cluster_properties.clustering_inds_agg{2},2);



[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('JDM25','JDM25caff',mocapmasterdirectory);
[mocapstruct_caff2] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);
ML_features_caff2 = get_supervised_features(mocapstruct_caff2,mocapstruct_caff2.modular_cluster_properties.clustering_inds_agg{2},2);

subset = 1:100:100*20000;

 plot_feature_spaces(ML_features,subset,'b')
  plot_feature_spaces(ML_features_prel,subset,'g')
plot_feature_spaces(ML_features_amph,subset,'r')
  plot_feature_spaces(ML_features_lesion,subset,'r')

  plot_feature_spaces(ML_features_caff2,subset,'r')

figure(55)
plot3(ML_features.pose_score(subset,1),ML_features.pose_score(subset,2),ML_features.pose_score(subset,4),'+b','MarkerSize',1)
hold on
plot3(ML_features_caff2.pose_score(subset,1),ML_features_caff2.pose_score(subset,2),ML_features_caff2.pose_score(subset,4),'+r','MarkerSize',1)
hold off

[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('JDM25','JDM25bipost',mocapmasterdirectory);
[mocapstruct_concatenated] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);


%% jdm 32
[descriptor_struct_6,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32postbi',mocapmasterdirectory);
[mocapstruct_post_bi] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_6,mocapfiletimes,0,0);
mocapstruct_post_bi = mocapstruct_post_bi.mocap_struct;




% 
% subset = 1:100:100*20000;%size(ML_features.pose_score,1);
% figure(55)
% plot3(ML_features.pose_score(subset,1),ML_features.pose_score(subset,2),log(ML_features.dyn_score(subset,2)),'+','MarkerSize',1)
% 
% figure(535)
% plot3(ML_features.appearance_features(subset,1),ML_features.appearance_features(subset,2),ML_features.appearance_features(subset,3),'+','MarkerSize',1)
% 
% plot3(ML_features.appearance_features_agg_score_whitened(subset,4),...
%     ML_features.appearance_features_agg_score_whitened(subset,2),...
%     ML_features.appearance_features_agg_score_whitened(subset,3),'+','MarkerSize',1)
% 
% figure(539)
% plot3(ML_features.dyn_score(subset,1),ML_features.dyn_score(subset,2),ML_features.dyn_score(subset,3),'+','MarkerSize',1)
% plot3(log(ML_features.dyn_score(subset,2)),log(ML_features.dyn_score(subset,3)),log(ML_features.dyn_score(subset,4)),'+','MarkerSize',1)
% 
% %mappedX = tsne(ML_features.pose_score(subset,:));
% mappedX = tsne(cat(2,ML_features.pose_score(subset,1:8),ML_features.appearance_features_agg_score_whitened(subset,1:6)));
% 
% 
% 

% figure(57)
% plot3(mappedX(:,1),mappedX(:,2), ML_features.pose_window_sd(1,subset,1),'+')
% ML_features.pose_score

% 
% figure(58)
% dynamics_inds = 1:100:100*23000;
% plot3(ML_features.spectrogram_pcs_trunk(dynamics_inds,1),ML_features.spectrogram_pcs_trunk(dynamics_inds,2),...
%     ML_features.spectrogram_pcs_trunk(dynamics_inds,3),'+','MarkerSize',1)
% 
% figure(588)
% dynamics_inds = 1:100:100*23000;
% plot3(ML_features.spectrogram_pcs_head(dynamics_inds,1),ML_features.spectrogram_pcs_head(dynamics_inds,2),...
%     ML_features.spectrogram_pcs_head(dynamics_inds,3),'+','MarkerSize',1)
% 
% 

subset = 1:100:100*5000;%size(ML_features.pose_score,1);
mappedX = tsne(cat(2,ML_features.pose_score(subset,1:6),ML_features.appearance_features_agg_score_whitened(subset,1:6)));


mappedX_dyn_angle =  tsne(cat(2,ML_features.spectrogram_pcs_head(subset,1:15),ML_features.spectrogram_pcs_trunk(subset,1:15)));
mappedX_dyn_head=  tsne(cat(2,ML_features.spectrogram_pcs_head_angle(subset,1:15)));


mappedX_dyn =  tsne(cat(2,ML_features.spectrogram_pcs_head(subset,1:5),ML_features.spectrogram_pcs_trunk(subset,1:5),...
   ML_features.spectrogram_pcs_hipL(subset,1:5),ML_features.spectrogram_pcs_hipR(subset,1:5)));

mappedX_joint =  tsne(cat(2,ML_features.pose_score(subset,1:10),ML_features.appearance_features_agg_score_whitened(subset,1:6),...
    ML_features.spectrogram_pcs_trunk_angle(subset,1:15),ML_features.spectrogram_pcs_head_angle(subset,1:15)));


figure(665)
plot(mappedX(:,1),mappedX(:,2),'+')


figure(666)
plot(mappedX_dyn_angle(:,1),mappedX_dyn_angle(:,2),'+')


figure(666)
plot(mappedX_dyn(:,1),mappedX_dyn(:,2),'+')

figure(667)
plot(mappedX_joint(:,1),mappedX_joint(:,2),'+')


framesubset = intersect(mocapstruct_caff.modular_cluster_properties.clustering_inds_agg{2},mocapstruct_caff.modular_cluster_properties.clipped_index{2});
framesubset_true = intersect(mocapstruct_caff.modular_cluster_properties.clustering_inds_agg{2},mocapstruct_caff.modular_cluster_properties.clipped_index{2});

framesubset = framesubset(subset);

% 
% 
% middleposes = framesubset(intersect((intersect(find(mappedX(:,1)>-50),find(mappedX(:,2)>50))),...
%     (intersect(find(mappedX(:,1)<300),find(mappedX(:,2)>103)))));
% 
% subsethere = framesubset(intersect(find(mappedX(:,1)>-200),find(mappedX(:,2)<-100)));
% animate_markers_aligned_fullmovie(mocapstruct_caff,subsethere)

%filename_mc = mocapfilearray{1};

%write commands here to concatenate the filearrays etc. together

%[mocapstruct_post] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes);

mocapstruct = mocapstruct_concatenated;%.mocap_struct;

%compare the eigen clusters based on specific timepoins
cluster_here = [8];
[modular_cluster_properties] = get_modularclustproperties(mocapstruct);
mocapstruct.rest_frames = 1:10;
mocapstruct.move_frames = 11:3000000;
%% this gives you the frames when particular m
[modular_cluster_properties] = get_clustering_features(mocapstruct,modular_cluster_properties,cluster_here)  ;

%% to try: SPINEF, SPINEL , HEAD(aligned), Or full pose (modular set 2, inc. the hips and knees and offset)

%visualize
animate_markers_aligned_fullmovie(mocapstruct,modular_cluster_properties.clustering_inds_agg{2}(1:20:end))

% start te gui
cd 'C:\Users\Jesse Marshall\Documents\GitHub\Movement_analysis\FATGUI-master'

 
 
 
MCC;

%% load the annotation structures
outputstruct_bipost = load('OutputStruct_bilDLS.mat');
outputstruct_bipost2 = load('OutputStruct_bidls_2.mat');
outputstruct_mcpost = load('OutputStruct.mat');

outputstruct_caff = load('OutputStruct_Vicon8_caff.mat');

fulloutput_DLS_cat = mergeStructs_JDM(outputstruct_bipost.output.GlobalBehavior,outputstruct_bipost2.output.GlobalBehavior);



%% select the behavior to examine
compare_string = 'LBodyGroom';
behavior_list = {'Walk','WetDogShake','FaceWipe','RearUp'};
behavior_list = {'RHeadScratch','LArmScratch','RArmScratch','Anogenitalgroom','RBodyGroom'};
behavior_list = {'Sniffstill'};

mocaplist = {mocapstruct_concatenated,mocapstruct_caff,mocapstruct_post_bi};   
fulloutput_DLS_cat = mergeStructs_JDM(outputstruct_bipost.output.Head,outputstruct_bipost2.output.Head);

%outputlist = {fulloutput_DLS_cat,outputstruct_caff.output.GlobalBehavior,outputstruct_mcpost.output.GlobalBehavior};
outputlist = {fulloutput_DLS_cat,outputstruct_caff.output.Head,outputstruct_mcpost.output.Head};

mocaptag = {'bi dls','caff','bi mc'};

savedirectory_here = strcat(savedirectory,'Annotationvideotests',filesep);
mkdir(savedirectory_here);
colorarray = {'r','b','g','k','c','y'};
dovideo = 0;
for mm = 1:4%numel(behavior_list)

    for kk = [1,3]%1:numel(mocaplist)
        if ~isempty(outputlist{kk}.(behavior_list{mm}))
[fulloutput,indivbouts] = fillannotationgaps(outputlist{kk}.(behavior_list{mm}),   20);
if (dovideo)
v = VideoWriter(strcat(savedirectory_here,behavior_list{mm},'_',mocaptag{kk}),'MPEG-4');
                open(v) 
                outputlength = min(numel(fulloutput),9000);
M = animate_markers_aligned_fullmovie(mocaplist{kk},fulloutput(1:10:outputlength));
 writeVideo(v,M)
 close(v)
end
if numel(fulloutput)>300
 multi_plot_marker_characteristics_timerange(mocaplist{kk},fulloutput,colorarray{kk},kk)
end
        end
        
if (kk == numel(mocaplist))
figure(200)
legend(mocaptag)
print('-dpng',strcat(savedirectory_here,'psdplot_',behavior_list{mm},'.png'))
close figure 200
 end
    
    end
end
 


compare_string = 'FaceWipes';

subset_of_frames_annotated = 1:100000;

[fulloutput_caff,indivbouts_caff] = fillannotationgaps(outputstruct_caff.output.GlobalBehavior.(compare_string),   20);
[fulloutput_caff_struct,indivbouts_caff_struct] = fillannotationgaps_struct(outputstruct_caff.output.GlobalBehavior,20);
% get annotated frames and frames where we have features
[outputvector,fieldnames_beh,observed_behaviors] = generate_categorical_output(fulloutput_caff_struct,numel(subset_of_frames_annotated));
fieldnames_beh = cat(1,'null',fieldnames_beh);

%make_ethogram(annotation_vec,fieldnames_beh)

annotated_subset_full = outputvector(intersect(framesubset_true,subset_of_frames_annotated));
[~,annotated_subset_good] = intersect(framesubset_true,subset_of_frames_annotated);

[candidate_frames,score] = findsimilarframes_mlfeatures_multifeatures(ML_features,annotated_subset_full,annotated_subset_good,1:2*10^6);

%[outputvector2,fieldnames_beh,observed_behaviors] = generate_categorical_output(candidate_frames,numel(subset_of_frames_annotated));

unique_predictions = unique(candidate_frames(~isnan(candidate_frames)));
fieldnames_beh(unique_predictions+1)

candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh ,'WetDogShake')==1)-1));
predicted_annotatedframes= candidate_frames_vec(candidate_frames_vec<80000);
predicted_unannotatedframes= candidate_frames_vec(candidate_frames_vec>80000);

animate_markers_aligned_fullmovie(mocapstruct_caff,framesubset_true(predicted_unannotatedframes(1:10:end))')
animate_markers_aligned_fullmovie(mocapstruct_caff,framesubset_true(predicted_annotatedframes(1:10:end))')

%% plot on the tsne
colors_plot = hsv(80);
colors_plot = colorcube(80)./2+lines(80)./2;
legendnames = cell(1,0);
for mm =1:max(candidate_frames)
   [ framesubset_cand,framesubset_tsne,~] = intersect(subset,find(candidate_frames == mm));
   if numel(framesubset_cand)
figure(559)
plot3(ML_features.pose_score(framesubset_cand,1),ML_features.pose_score(framesubset_cand,2),...
    ML_features.pose_score(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
hold on

% figure(56900)
% plot(mappedX_dyn(framesubset_tsne,1),mappedX_dyn(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
% hold on
 legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};
% 
% 
% figure(569)
% plot(mappedX(framesubset_tsne,1),mappedX(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
% hold on


figure(570)
plot(mappedX_joint(framesubset_tsne,1),mappedX_joint(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
hold on

% 
% 
% figure(572)
% plot(mappedX_dyn_angle(framesubset_tsne,1),mappedX_dyn_angle(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
% hold on




figure(571)
plot3(ML_features.spectrogram_pcs_head(framesubset_cand,1),ML_features.spectrogram_pcs_head(framesubset_cand,2),...
    ML_features.spectrogram_pcs_head(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
hold on

   end
end
figure(56900)
legend(legendnames)
    figure(559)
legend(legendnames)
    figure(569)
legend(legendnames)
    figure(570)
legend(legendnames)
   figure(571)
legend(legendnames)
  figure(572)
legend(legendnames)

make_ethogram(outputvector,fieldnames_beh)
 %   make_dotplot(outputvector,fieldnames_beh,ML_features)
        make_dotplot(    candidate_frames,fieldnames_beh,ML_features)

    %make_transitionplot
    %get frequency of each behavior
    

    %% look at tsne of a subset
    candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh , 'LBodyGroom')==1)-1));
%multi_plot_marker_characteristics_timerange(mocapstruct_caff,framesubset_true(candidate_frames_vec),'c',1)
[fulloutput,indivbouts] = fillannotationgaps(framesubset_true(candidate_frames_vec),20);
indivbouts_num = cellfun(@numel,indivbouts);
indivbouts(indivbouts_num<30) = [];
plot_behavior_examples(mocapstruct_caff,fulloutput,indivbouts);


    frames_use = candidate_frames_vec(1:10:end);
    frames_use = frames_use(1:min(numel(frames_use),10000));
mappedX_sub = tsne(ML_features.pose_score(frames_use,:));
figure(46)
plot(mappedX_sub(:,1),mappedX_sub(:,2),'+')





middleposes = (intersect((intersect(find(mappedX_sub(:,1)>-20),find(mappedX_sub(:,2)<20))),...
    (intersect(find(mappedX_sub(:,1)<0),find(mappedX_sub(:,2)>0)))));

animate_markers_aligned_fullmovie(mocapstruct_caff,framesubset_true((frames_use(middleposes)))')

    
%% do for a single feature
[~,fullframescaff_good] = intersect(framesubset_true,fulloutput_caff);

[candidate_frames,score] = findsimilarframes_mlfeatures(ML_features,fullframescaff_good,annotated_subset_good,1:10^6);

%[~,candidateframes] = intersect(framesubset_true,1:2000000);
candidateframes = find(candidate_frames == 1);
candidateframes_score = find(score(:,2)>0.5);
predicted_annotatedframes = candidateframes(candidateframes<80000);

predicted_unannotatedframes= candidateframes(candidateframes>80000);
predicted_unannotatedframes_sc = candidateframes_score(candidateframes_score>80000);


animate_markers_aligned_fullmovie(mocapstruct_caff,framesubset_true(predicted_unannotatedframes(1:10:end))')
animate_markers_aligned_fullmovie(mocapstruct_caff,framesubset_true(fullframescaff_good(1:10:end))')
animate_markers_aligned_fullmovie(mocapstruct_caff,framesubset_true(predicted_annotatedframes(1:10:end))')
animate_markers_aligned_fullmovie(mocapstruct_caff,framesubset_true(predicted_unannotatedframes_sc(1:10:end))')


%animate_markers_aligned_fullmovie(mocapstruct_caff,fulloutput_caff(1:10:end)')

[~,subset_1] = intersect(framesubset_true,fulloutput_caff);
[subset_2] = (predicted_unannotatedframes);

figure(55)
plot3(ML_features.pose_score(subset,1),ML_features.pose_score(subset,2),ML_features.pose_score(subset,3),'+b','MarkerSize',1)
hold on
plot3(ML_features.pose_score(subset_2,1),ML_features.pose_score(subset_2,2),ML_features.pose_score(subset_2,3),'+r','MarkerSize',1)
plot3(ML_features.pose_score(subset_1,1),ML_features.pose_score(subset_1,2),ML_features.pose_score(subset_1,3),'+g','MarkerSize',1)
hold off

[~,tsnesubset] = intersect(subset,subset_1);
[~,tsnesubset_rf] = intersect(subset,subset_2);

figure(56)
plot(mappedX(:,1),mappedX(:,2),'b+')
hold on
plot(mappedX(tsnesubset,1),mappedX(tsnesubset,2),'r+')
plot(mappedX(tsnesubset_rf,1),mappedX(tsnesubset_rf,2),'g+')

hold off


[fulloutput_DLS2,indivbouts_dls] = fillannotationgaps(fulloutput_DLS_cat.(compare_string),20);
[fulloutput_mc,indivbouts_mc] = fillannotationgaps(outputstruct_mcpost.output.GlobalBehavior.(compare_string),   20);

 
%animate_markers_aligned_fullmovie(mocapstruct_concatenated,fulloutput_DLS(1:10:end))
plot_behavior_examples(mocapstruct_caff,fulloutput_caff,indivbouts_caff);


%% plot the marker characteristics during this behavior
legend_tags = {'caff','bidls','bimc'};
cellstruct = cell(1,3);
cellstruct{1} = mocapstruct_caff;
cellstruct{2} = mocapstruct_concatenated;
cellstruct{3} = mocapstruct_post_bi;
timestruct{1} = fulloutput_caff;
timestruct{2} = fulloutput_DLS2;
timestruct{3} = fulloutput_mc;


for kk =1:3
multi_plot_marker_characteristics_timerange(cellstruct{kk},timestruct{kk},colorarray{kk},kk)
end
% figure(200)
% subplot(2,1,1)
% legend(legend_tags)
% subplot(2,1,2)
% legend(legend_tags)
% % save here
% figure(306)
% legend(legend_tags)
% figure(307)
% subplot(1,2,1)
% legend(legend_tags)
% subplot(1,2,2)
% legend(legend_tags)



%compare_plot_marker_characteristics_timerange(mocapstruct_post_bi,fulloutput_mc,mocapstruct_caff,fulloutput_caff)

%candidateframes = findsimilarframes(mocapstruct_concatenated,fulloutput_DLS2,[1,2,3,4,7,8,13,14],1:500000);
%animate_markers_aligned_fullmovie(mocapstruct_concatenated,candidateframes(1:10:end)')



