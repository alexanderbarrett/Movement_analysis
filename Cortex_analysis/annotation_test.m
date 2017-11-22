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

%% write commands here to concatenate the filearrays etc. together

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

testframes = 100000:200000;

outputstruct_bipost = load('OutputStruct_bilDLS.mat');
outputstruct_caff = load('OutputStruct_Vicon8_caff.mat');

fulloutput_caff = fillannotationgaps(outputstruct_caff.output.GlobalBehavior.FaceWipes,   20);
fulloutput_DLS = fillannotationgaps(outputstruct_bipost.output.GlobalBehavior.Walk,   20);

% 
% output2 = load( 'Backup_OutPutStruct\Backup_OutputStruct_November-16-2017-14-21-33.mat')
% fulloutput = fillannotationgaps(outputstruct.output.GlobalBehavior.Walk,20);
% fulloutput2 = fillannotationgaps(outputstruct.output.GlobalBehavior.Walk,20);
% 
% fulloutput= fillannotationgaps(outputstruct.output.GlobalBehavior.RearUp,10);
% fulloutput= fillannotationgaps(outputstruct.output.GlobalBehavior.WetDogShake,   20);
% 
% fulloutput= fillannotationgaps(outputstruct.output.GlobalBehavior.RBodyGroom,   20);
% fulloutput= fillannotationgaps(outputstruct.output.GlobalBehavior.FaceWipes,   20);

animate_markers_aligned_fullmovie(mocapstruct_caff,fulloutput_caff(1:10:end))
animate_markers_aligned_fullmovie(mocapstruct_concatenated,fulloutput_DLS(1:10:end))


figure(77)
subplot(2,1,1)
plot(mocapstruct_caff.markers_aligned_preproc.HeadF(fulloutput_caff,3),'b')
hold on
plot(mocapstruct_caff.markers_aligned_preproc.ShoulderR(fulloutput_caff,3),'r')
hold off

subplot(2,1,2)
plot(mocapstruct_concatenated.markers_aligned_preproc.HipL(fulloutput_DLS,3),'b')
hold on
plot(mocapstruct_concatenated.markers_aligned_preproc.HipR(fulloutput_DLS,3),'r')
hold off



plot(mocapstruct.markers_aligned_preproc.ElbowR(fulloutput,3),'c')
plot(mocapstruct.markers_aligned_preproc.ElbowL(fulloutput,3),'g')

plot(mocapstruct.markers_aligned_preproc.HeadF(fulloutput,3),'r')
%plot(mocapstruct.markers_aligned_preproc.Offset1(fulloutput,3),'c')



