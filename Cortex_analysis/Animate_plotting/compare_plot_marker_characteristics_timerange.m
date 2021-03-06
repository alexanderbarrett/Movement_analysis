function compare_plot_marker_characteristics_timerange(mocapstruct,timerange,mocapstruct2,timerange2)

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
marker_velocity_agg2 = cell(1,numel(mocapstruct.markernames));
xval_agg = cell(1,numel(mocapstruct.markernames));

for ll = 1:numel(mocapstruct.markernames)

    marker_here = struct('singlemarker',[]);

[~,badframeintersect,~]= intersect(timerange,mocapstruct.bad_frames_agg{ll});
goodframes = timerange(setxor(1:numel(timerange),badframeintersect));
marker_here.singlemarker = mocapstruct.markers_preproc.(mocapstruct.markernames{ll});
params.fps = fps;
[marker_clipped,clipped_index] = hipass_clip_fragments(marker_here,goodframes,params);
% figure(44)
% plot(marker_clipped.singlemarker)
pwelch_no = 1;


P_sum = [];
for mm = 1:3
[Pxx,fxx] = pwelch((marker_clipped.singlemarker(:,mm)),fps*pwelch_no,floor(0.5*fps),fps*pwelch_no,fps,'onesided');
if (mm == 1)
    P_sum = Pxx;
else
   P_sum = P_sum + Pxx; 
end
end
% 
% marker_here2.singlemarker = mocapstruct2.markers_preproc.(mocapstruct2.markernames{ll})(timerange2,:);
% params.fps = fps;
% [~,badframeintersect2,~]= intersect(timerange2,mocapstruct.bad_frames_agg{ll});
% 
     marker_here2 = struct('singlemarker',[]);
[~,badframeintersect,~]= intersect(timerange2,mocapstruct2.bad_frames_agg{ll});
goodframes2 = timerange2(setxor(1:numel(timerange2),badframeintersect));
marker_here2.singlemarker = mocapstruct2.markers_preproc.(mocapstruct2.markernames{ll});
params.fps = fps;
[marker_clipped2,clipped_index] = hipass_clip_fragments(marker_here2,goodframes2,params);
 figure(44)
plot(marker_clipped2.singlemarker)

%[marker_clipped2,clipped_index] = hipass_clip(marker_here2,badframeintersect2,params);
P_sum2 = [];
for mm = 1:3
[Pxx2,fxx] = pwelch((marker_clipped2.singlemarker(:,mm)),fps*pwelch_no,floor(0.5*fps),fps*pwelch_no,fps,'onesided');
if (mm == 1)
    P_sum2 = Pxx2;
else
   P_sum2 = P_sum2 + Pxx2; 
end
end


figure(200)
subplot_tight(subplot_row,subplot_cols,ll)
plot(fxx,log10(P_sum),'b')
hold on
plot(fxx,log10(P_sum2),'r')
hold off
xlim([1 50])
ylabel('SP ')
xlabel('Freq')
title(mocapstruct.markernames{ll})
if ll == 4
legend('input1','input2','location','North')
end
print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerPSDs_move_tr.png'))

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
% first dataset
veltemp =  diff(marker_clipped.singlemarker(:,1),6).^2;
for jj = 2:3
    veltemp = veltemp+diff(marker_clipped.singlemarker(:,jj),6).^2;
end
marker_velocity = sqrt(veltemp./3);
marker_velocity_agg{ll} = marker_velocity;

%% second dataset
veltemp2 =  diff(marker_clipped2.singlemarker(:,1),6).^2;
for jj = 2:3
    veltemp2 = veltemp2+diff(marker_clipped2.singlemarker(:,jj),6).^2;
end
marker_velocity2 = sqrt(veltemp2./3);
marker_velocity_agg2{ll} = marker_velocity2;


xval_agg{ll} = ll.*ones(1,numel(marker_velocity));

%% plot here
  figure(300)
 subplot_tight(subplot_row,subplot_cols,ll)

    [n,x] = hist( marker_velocity,0:0.1:10);
  %  bar(x,log10(n),'b','EdgeColor','none');
       semilogy(x,log10(n),'b','Linewidth',2);

     hold on
     [n2,x2] = hist( marker_velocity2,0:0.1:10);
   % bar(x2,log10(n2),'r','EdgeColor','none');
        semilogy(x2,log10(n2),'r','Linewidth',2);

    hold off
%     [n2,x2] = hist(marker_velocity,0:0.1:5);
%     bar(x2,log10(n2),'r');
%     hold off
 title(mocapstruct.markernames{ll})
    xlim([0 10])
    xlabel('Marker Velociy [mm/frame]')
    ylabel('Fraction of bins')
  print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerVelocities_move_tr.png'))


  
end

figure(100)
bar(30/6.*cellfun(@mean,marker_velocity_agg))
box off
set(gca,'XTick',1:numel(mocapstruct.markernames),'XTickLabels',mocapstruct.markernames)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Marker Velocity (cm/s)')
xlim([0 numel(mocapstruct.markernames)+1])


figure(101)
bar(30/6.*cellfun(@mean,marker_velocity_agg2))
box off
set(gca,'XTick',1:numel(mocapstruct2.markernames),'XTickLabels',mocapstruct2.markernames)
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
