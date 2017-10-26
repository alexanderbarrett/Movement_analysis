%examine hand cluster quality

figure(11)
subplot(2,1,1)
plot(mocapstruct.markers_aligned_preproc.ArmL)
hold on
frames_tot = zeros(1,size(mocapstruct.markers_aligned_preproc.ArmL,1));
frames_tot(mocapstruct.bad_frames_agg{12}) = 1;
plot(50*frames_tot,'g')
hold off
subplot(2,1,2)
frame_axis 
plot(mocapstruct.markers_aligned_preproc.ArmL(find(frames_tot==0),:))
