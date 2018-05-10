function   plot_annotation_params(mocapstruct,annotation_struct,fieldname)



fulloutput_DLS2 = fillannotationgaps(fulloutput_DLS_cat.FaceWipe,20);




params.fps = 300;
markerout2 = hipass_clip_fragments(mocapstruct_caff.markers_aligned_preproc ,fulloutput_caff,params);
xlimhere = 2000;
figure(77)
subplot(2,1,1)
plot(markerout2.ElbowR(:,3),'b')
hold on
plot(markerout2.HeadF(:,3),'r')
hold off
xlim([0 xlimhere])

markerout = hipass_clip_fragments(mocapstruct_concatenated.markers_aligned_preproc ,fulloutput_DLS2,params);

subplot(2,1,2)
plot(markerout.ShoulderR(:,3),'b')
hold on
plot(markerout.HeadF(:,3),'r')
hold off
xlim([0 xlimhere])

end
