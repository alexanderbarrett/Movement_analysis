%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots_videoexamples',filesep);
mkdir(savedirectory);



%% load or create struct
createmocapfilestruct('Vicon3',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct = loadmocapfilestruct('Vicon3',mocapmasterdirectory);
mocapfilestruct.mocapdir = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\';
mocapfilestruct.mocapdir = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\';

[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon3','Vicon3_htr',mocapmasterdirectory);
mocapfilestruct.mocapdir = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\';


filearray = {'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid1.htr',...
 'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid2.htr',...
    'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid3.htr',...
    'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid4.htr',...
    'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid5.htr',...
    };

    [htr_struct] = concat_htr(filearray);
    segnames = fieldnames(htr_struct.tree);

    rotations_quat_exp = cell(1,numel(mocapstruct_pre.rotations));

for mm = limbstart:limbstop
    initial_sign = htr_struct.agg_rotations_quat{mm}(1,:);
      figure(440+mm)
 ax1=   subplot(3,1,1)
plot(htr_struct.agg_rotations_quat{mm})
    title(segnames{mm})
  
    sign_product = sum(bsxfun(@times,(htr_struct.agg_rotations_quat{mm}(:,2:4)),initial_sign(2:4)),2);
    ax2=   subplot(3,1,2)
plot(sign_product)

  ax3=   subplot(3,1,3)
plot(bsxfun(@times,htr_struct.agg_rotations_quat{mm},sign(sign_product)))
    title(segnames{mm})
%bsxfun(@times,htr_struct.agg_rotations_quat{mm},[1 -1 -1 -1]))
    title(segnames{mm})
%   ax2=  subplot(3,1,2)
%        % rotations_quat_exp{mm} = bsxfun(@times,htr_struct.agg_rotations_quat{mm}(:,1:3),...
%        %     (2*acos(htr_struct.agg_rotations_quat{mm}(:,4))./(sqrt(sum(htr_struct.agg_rotations_quat{mm}(:,1:3).^2,2) ))));
%  test=SpinConv('QtoEA321',htr_struct.agg_rotations_quat{mm}(:,:));
%  plot(test)
%  
%    ax3=  subplot(3,1,3)
%  test=SpinConv('EA321toQ',test);
%  plot(test)

%plot(rotations_quat_exp{mm})
linkaxes([ax1,ax2,ax3],'x')

end    

figure(444)
subplot(2,1,1)
plot(quatnorm(htr_struct.agg_rotations_quat{7}))
hold on
plot(quatmod(htr_struct.agg_rotations_quat{7}))
subplot(2,1,2)
plot(quatconj((htr_struct.agg_rotations_quat{7})))

%% code for looking at joint angles
[mocapstruct_pre] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,0,1,filearray);

good_joints = setxor(1:10,[3,6,8]);
corresponding_markers = [1:10,17,18];
bad_frames_joints = unique(cat(2,mocapstruct_pre.bad_frames_agg{corresponding_markers}));
good_frames = setxor(1:size(mocapstruct_pre.agg_rotations_quat{1},1),bad_frames_joints);

rot_quat_exp = cell(1,18);
agg_features = [];
for mm = good_joints
    
    plot_here = mocapstruct_pre.agg_rotations_quat{mm};
    plot_here(bad_frames_joints,:) = nan;
figure(440+mm)
ax1 = subplot(3,1,1)
plot(mocapstruct_pre.agg_rotations_quat{mm})

%% show the frames without missing pts
ax2 = subplot(3,1,2)

plot( plot_here)
title(mocapstruct_pre.jointnames{mm})

%% exponentially map the quaternions to get a representation in R3
 rotations_quat_exp{mm} = bsxfun(@times,mocapstruct_pre.agg_rotations_quat{mm}(:,2:4),...
             (2*acos(mocapstruct_pre.agg_rotations_quat{mm}(:,1))./(sqrt(sum(mocapstruct_pre.agg_rotations_quat{mm}(:,2:4).^2,2) ))));
ax3 = subplot(3,1,3)
plot( rotations_quat_exp{mm})
title(mocapstruct_pre.jointnames{mm})

linkaxes([ax1,ax2,ax3],'x')

%% these are the features to use -- median filter and remove means
figure(550+mm)
clipped_joint_angle = struct('singlemarker',rotations_quat_exp{mm});
params.fps = 300;
clipped_joint_angle_out = hipass_clip_fragments(clipped_joint_angle,good_frames,params);
plot(medfilt1(clipped_joint_angle_out.singlemarker,5))
title(mocapstruct_pre.jointnames{mm})

%% create aggregate features 
 agg_features = cat(1,agg_features, medfilt1(clipped_joint_angle_out.singlemarker,5)',agg_features);
% opts.fps =300;
% opts.clustering_window = opts.fps./2;
% opts.clustering_overlap = opts.fps./4;
% [dyad_out,fr,time_clustering] = get_dyadic_spectrogram(agg_features,opts);
% imagesc(dyad_out)

end

figure(44)
plot(agg_features')




%% look at the tait-bryant representation to see yaw pitch and roll of the joints
limbhere = 5;
figure(5556)

 title(segnames{limbhere})
legend('yaw','pitch','roll')



rotations_eul = cell(1,numel(mocapstruct_pre.rotations));
rotations_eul_interp = cell(1,numel(mocapstruct_pre.rotations));

rotations_quat = cell(1,numel(mocapstruct_pre.rotations));
rotations_quat_exp = cell(1,numel(mocapstruct_pre.rotations));
rotations_quat_interp = cell(1,numel(mocapstruct_pre.rotations));


for mm = limbstart:limbstop

    %do interpolated
    rotations_eul_interp{mm} = rotm2eul(mocapstruct_pre.rotations{mm});
[nanindsall] = find(isnan((rotations_eul_interp{mm})));
[naninds,~] = ind2sub(size(rotations_eul_interp{mm}),nanindsall);
naninds = unique(naninds);
%rotations{mm}(isnan(rotations{mm})) = 0;
rotations_quat_interp{mm} = rotm2quat(mocapstruct_pre.rotations{mm});
rotations_quat_interp{mm}(naninds,:) = nan;
       rotations_quat_interp{mm}(find( sum(sum(mocapstruct_pre.rotations{mm},1),2)==0),:) = nan;

    
    
    rotations_eul{mm} = rotm2eul(rotations{mm});
[nanindsall] = find(isnan((rotations_eul{mm})));
[naninds,~] = ind2sub(size(rotations_eul{mm}),nanindsall);
naninds = unique(naninds);
rotations{mm}(isnan(rotations{mm})) = 0;
rotations_quat{mm} = rotm2quat(rotations{mm});
rotations_quat{mm}(naninds,:) = nan;

    rotations_quat_exp{mm} = bsxfun(@times,rotations_quat{mm}(:,1:3),(2*acos(rotations_quat{mm}(:,4))./(sqrt(sum(rotations_quat{mm}(:,1:3).^2,2) ))));
rotations_quat_exp{mm}(naninds,:) = nan;
   rotations_quat_exp{mm}(find( sum(sum(rotations{mm},1),2)==0),:) = nan;
   
    figure(44+mm)
    title(segnames{mm})
  ax1 =  subplot(3,1,1);
  
plot(rotations_quat{mm})
 ax2 = subplot(3,1,2);
plot(rotations_quat_exp{mm})
 ax3 = subplot(3,1,3);
plot(rotations_quat_interp{mm})
linkaxes([ax1,ax2,ax3],'x')


end

for mm = limbstart:limbstop
end

figure(444)
plot(rotations_quat_exp{1})



[modular_cluster_properties_social] = get_modularclustproperties(mocapstruct_social);