%% annotation test

%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
%mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
loadmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server

mocapfilestruct = loadmocapfilestruct('JDM25',mocapmasterdirectory);

%alternate: Vicon8_caff, Vicon8_dlslesion
[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);
[mocapstruct_caff] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,1,1);

[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('JDM25','JDM25bipost',mocapmasterdirectory);
[mocapstruct_concatenated] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);


%% jdm 32
[descriptor_struct_6,mocapfilearray,mocapfilestruct_JDM32,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('JDM32','JDM32postbi',mocapmasterdirectory);
[mocapstruct_post_bi] = preprocess_mocap_data(mocapfilearray,mocapfilestruct_JDM32,descriptor_struct_6,mocapfiletimes,0,0);
mocapstruct_post_bi = mocapstruct_post_bi.mocap_struct;
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
compare_string = 'Walk';
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
 
[fulloutput_caff,indivbouts_caff] = fillannotationgaps(outputstruct_caff.output.GlobalBehavior.(compare_string),   20);

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



