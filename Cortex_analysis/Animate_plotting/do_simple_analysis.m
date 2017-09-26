function do_simple_analysis(traceuse,veluse,region_plot,tag,fps,markers_preproc,plotdirectory,markercolor,links)


  %  traceuse =traceuse(region_plot);
    %    traceuse =traceuse(region_plot);

% figure(33)
% %plot(0:1./300:((size(veluse,2)-1)./300),squeeze(marker_velocity(markeruse, :,4)),'b')
% plot(0:1./fps:((size(traceuse,2)-1)./fps),0.1*squeeze((traceuse)),'b')
% 
% hold on
% temp = zeros(size(veluse));
% temp(bad_frames_agg{markeruse}) = 1;
% 
% plot(0:1./fps:((size(traceuse,2)-1)./fps),temp.*10,'g')
% hold off
% ylabel('Head Vel (mm)')
% xlabel('Time (s)')


[Pxx,fxx] = pwelch(traceuse,fps*2,0.5,fps*2,fps,'onesided');
figure(34)
plot(fxx,log10(Pxx))
xlim([1 50])
ylabel('Log10 Spectral Power ')
xlabel('Frequency (Hz)')
title(tag)
print('-dpng',strcat(plotdirectory,tag,'psd.png'));

%xlim([3.910*10^6/300 3.950*10^6/300])


[~,f,t,spectrogram_here] = spectrogram(traceuse,2*fps,fps*0.5,fps*2,fps);

figure(39)
subplot(2,1,1)
plot(0:1./fps:((numel(traceuse)-1)./fps),traceuse);

xlim([0 (numel(traceuse)-1)./fps])

subplot(2,1,2)
imagesc(t,f,log10(spectrogram_here))
ylim([0 50])
caxis([-10 2])
print('-dpng',strcat(plotdirectory,'tag.png'));
xlim([0 (numel(traceuse)-1)./fps])

h=figure(40)
markernames = fieldnames(markers_preproc);
[markers_centered] = center_marker_struct(markers_preproc,fps);

M=[];
animate_markers_clusteraligned(markers_centered,markers_preproc,region_plot',markernames,markercolor,links,M,1,'1',h,tag);


% figure(40)
% marker_1 = 11;
% mean_1 = marker_velocity(5,:,4);%markers_preproc.(marker_names{ marker_1})(:,3);
% 
% %mean_1 = abs_velocity_antialiased(5,:);%markers_preproc.(marker_names{ marker_1})(:,3);
% mean_1(bad_frames_agg{ 5}) = nan;
% 
% plot_1 = marker_velocity(marker_1,:,4);%markers_preproc.(marker_names{ marker_1})(:,3);
% 
% plot_1 = abs_velocity_antialiased(marker_1,:);%markers_preproc.(marker_names{ marker_1})(:,3);
% plot_1(bad_frames_agg{ marker_1}) = nan;
% 
% marker_2 = 15;
% plot_2 = marker_velocity(marker_2,:,4);%markers_preproc.(marker_names{ marker_1})(:,3);
% 
% plot_2 = abs_velocity_antialiased(marker_2,:);%markers_preproc.(marker_names{marker_2 })(:,3);
% plot_2(bad_frames_agg{marker_2 }) = nan;
% 
% %subplot(2,1,2)
% %  plot(markers_preproc.(marker_names{1})(goodtimeshere,:))
% %plot(markers_preproc.(marker_names{ marker_1})(:,3),'c')
% plot(plot_1,'b')
% 
% hold on
% 
% plot(plot_2,'r')
% 
% %set(gca,'XTick',filestartpts(1,:),'XTickLabels',mocap_datecreate);
% %    filestartpts
% %mocap_datecreate
% hold off
% 
% nanmedian(plot_1)
% 
% nanmedian(plot_2)
% 
% badframes = cat(2,bad_frames_agg{ marker_1},bad_frames_agg{ marker_2});
% good_frames = setxor(1:numel(plot_1),badframes);
% 
% plot_1 = plot_1;
% plot_2 = plot_2;
% 
% figure(35)
% [C,LAGS] = xcorr(plot_1(good_frames)-nanmean(plot_1(good_frames)),plot_2(good_frames)-nanmean(plot_2(good_frames)),300,'Coeff');
% plot(LAGS./300,C)
% title(marker_names{markeruse})
% xlabel('time (s)')
% ylabel('acorr strength')
end



