function multi_plot_marker_characteristics_timerange(mocapstruct,timerange,colorin,iter)

fps= mocapstruct.fps;

% %% first plot the velocity histograms of each marker
% for mm = 1:size(marker_velocity,1)
%     length_vel_here = size(marker_velocity,2);
%   
% end
mocapstruct.markernames;

subplot_row = 4;

marker_velocity_agg = cell(1,numel(mocapstruct.markernames));
marker_velocity_agg2 = cell(1,numel(mocapstruct.markernames));
xval_agg = cell(1,numel(mocapstruct.markernames));


markersubset = [1,7,11,15,17,18];
subplot_cols = ceil(numel(markersubset)./subplot_row);

for ll = markersubset;%1:numel(mocapstruct.markernames)

    marker_here = struct('singlemarker',[]);

[~,badframeintersect,~]= intersect(timerange,mocapstruct.bad_frames_agg{ll});
goodframes = timerange(setxor(1:numel(timerange),badframeintersect));
marker_here.singlemarker = mocapstruct.markers_preproc.(mocapstruct.markernames{ll});
params.fps = fps;
[marker_clipped,clipped_index] = hipass_clip_fragments(marker_here,goodframes,params);
% figure(44)
% plot(marker_clipped.singlemarker)
pwelch_no = 1;


if numel(marker_clipped.singlemarker(:,1))>fps*pwelch_no
P_sum = [];
for mm = 1:3
[Pxx,fxx] = pwelch((marker_clipped.singlemarker(:,mm)),fps*pwelch_no,floor(0.5*fps),fps*pwelch_no,fps,'onesided');
if (mm == 1)
    P_sum = Pxx;
else
   P_sum = P_sum + Pxx; 
end
end


figure(200)
hold on
subplot_tight(subplot_row,subplot_cols,find(markersubset == ll))
plot(fxx,log10(P_sum),'color',colorin)
hold off
xlim([1 50])
ylabel('SP ')
xlabel('Freq')
ntitle(mocapstruct.markernames{ll},'location','North')
% if ll == 4
% legend('input1','input2','location','North')
% end
print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerPSDs_move_tr.png'))



%% plot the marker velocity
%marker_velocity = diff(marker_clipped.singlemarker(:,x),6).^2,jj);
% first dataset
veltemp =  diff(marker_clipped.singlemarker(:,1),6).^2;
for jj = 2:3
    veltemp = veltemp+diff(marker_clipped.singlemarker(:,jj),6).^2;
end
marker_velocity = sqrt(veltemp./3);
marker_velocity_agg{ll} = marker_velocity;



xval_agg{ll} = ll.*ones(1,numel(marker_velocity));

%% plot here
  figure(300)
subplot_tight(subplot_row,subplot_cols,find(markersubset == ll))
hold on
    [n,x] = hist( marker_velocity,0:0.1:10);
  %  bar(x,log10(n),'b','EdgeColor','none');
       semilogy(x,log10(n),'Color',colorin,'Linewidth',2);
set(gca,'yscale','log')
    hold off
%     [n2,x2] = hist(marker_velocity,0:0.1:5);
%     bar(x2,log10(n2),'r');
%     hold off
 ntitle(mocapstruct.markernames{ll},'location','north')
    xlim([0 10])
    xlabel('Marker Velociy [mm/frame]')
    ylabel('Fraction of bins')
  print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerVelocities_move_tr.png'))


end
end

figure(100+iter)
bar(30/6.*cellfun(@mean,marker_velocity_agg))
box off
set(gca,'XTick',1:numel(mocapstruct.markernames),'XTickLabels',mocapstruct.markernames)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Marker Velocity (cm/s)')
xlim([0 numel(mocapstruct.markernames)+1])




figure(120+iter)
boxplot(cat(1,marker_velocity_agg{:}),cat(2,xval_agg{:}),'Symbol','')
box off
set(gca,'XTick',1:numel(mocapstruct.markernames),'XTickLabels',mocapstruct.markernames)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Marker Velocity (cm/s)')
xlim([0 numel(mocapstruct.markernames)+1])
ylim([0 2])
  print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerVelocities_box_move.png'))
