function plot_marker_characteristics(mocapstruct)

fps= mocapstruct.fps;

% %% first plot the velocity histograms of each marker
% for mm = 1:size(marker_velocity,1)
%     length_vel_here = size(marker_velocity,2);
%   
% end

mocapstruct.markernames;

subplot_row = 4;
subplot_cols = ceil(numel(mocapstruct.markernames)./subplot_row);

marker_velocity_agg = cell(1,numel(mocapstruct.markernames));
xval_agg = cell(1,numel(mocapstruct.markernames));

for ll = 1:numel(mocapstruct.markernames)
marker_here = struct('singlemarker',[]);
marker_here.singlemarker = mocapstruct.markers_preproc.(mocapstruct.markernames{ll});



params.fps = fps;
[marker_clipped,clipped_index] = hipass_clip(marker_here,cat(2,mocapstruct.bad_frames_agg{ll},mocapstruct.rest_frames),params);



[Pxx,fxx] = pwelch(marker_clipped.singlemarker(:,3),fps*5,floor(0.5*fps),fps*5,fps,'onesided');


figure(200)
subplot_tight(subplot_row,subplot_cols,ll)
plot(fxx,log10(Pxx))
xlim([1 50])
ylabel('SP ')
xlabel('Freq')
title(mocapstruct.markernames{ll})
print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerPSDs_move.png'))

%xlim([3.910*10^6/300 3.950*10^6/300])


%[~,f,t,spectrogram_here] = spectrogram(marker_clipped.singlemarker(:,3),2*fps,fps*1,fps*2,fps);
% 
% figure(400+ll)
% %subplot(2,1,1)
% imagesc(t,f,log10(spectrogram_here))
% ylim([0 30])
% caxis([-10 0])
%

%% plot the marker velocity
%marker_velocity = diff(marker_clipped.singlemarker(:,x),6).^2,jj);
veltemp =  diff(marker_clipped.singlemarker(:,1),6).^2;
for jj = 2:3
    veltemp = veltemp+diff(marker_clipped.singlemarker(:,jj),6).^2;
end
marker_velocity = sqrt(veltemp./3);
marker_velocity_agg{ll} = marker_velocity;
xval_agg{ll} = ll.*ones(1,numel(marker_velocity));

  figure(300)
 subplot_tight(subplot_row,subplot_cols,ll)

    [n,x] = hist( marker_velocity,0:0.1:5);
    bar(x,log10(n),'b');
%     hold on
%     [n2,x2] = hist(marker_velocity,0:0.1:5);
%     bar(x2,log10(n2),'r');
%     hold off
 title(mocapstruct.markernames{ll})
    xlim([0 5])
    xlabel('Marker Velociy [mm/frame]')
    ylabel('Fraction of bins')
  print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerVelocities_move.png'))


end

figure(100)
bar(30/6.*cellfun(@mean,marker_velocity_agg))
box off
set(gca,'XTick',1:numel(mocapstruct.markernames),'XTickLabels',mocapstruct.markernames)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Marker Velocity (cm/s)')
xlim([0 numel(mocapstruct.markernames)+1])



figure(120)
boxplot(cat(1,marker_velocity_agg{:}),cat(2,xval_agg{:}),'Symbol','')
box off
set(gca,'XTick',1:numel(mocapstruct.markernames),'XTickLabels',mocapstruct.markernames)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Marker Velocity (cm/s)')
xlim([0 numel(mocapstruct.markernames)+1])
ylim([0 2])
  print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerVelocities_box_move.png'))
