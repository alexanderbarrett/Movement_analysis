function plotfractionmissing(mocapstruct)
marker_names = fieldnames(mocapstruct.markers_preproc);

fraction_missing_rest = zeros(1,numel(marker_names));
fraction_missing_move = zeros(1,numel(marker_names));

length_rest = numel(mocapstruct.rest_frames);
length_move = numel(mocapstruct.move_frames);

rest_dropout_lengths = cell(1,numel(marker_names));
move_dropout_lengths = cell(1,numel(marker_names));

fraction_missing_rest_nolong = zeros(1,numel(marker_names));
fraction_missing_move_nolong = zeros(1,numel(marker_names));

long_thresh = 15;

%size_mocap = size( mocapstruct.markers_preproc.(marker_names{1}),1);
for ll = 1:numel(marker_names)
   fraction_missing_rest(ll) = numel(find(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.rest_frames,:),2) == 0))./...
       length_rest; 
   fraction_missing_move(ll) = numel(find(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.move_frames,:),2) == 0))./...
       length_move; 
   
   %% get the distribution of missing times
   trace_to_label = zeros(1,size(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.rest_frames,:),1));
   trace_to_label(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.rest_frames,:),2)==0) = 1;   
   pixellist = bwconncomp(trace_to_label);
rest_dropout_lengths{ll} = cellfun(@numel,pixellist.PixelIdxList);


   trace_to_label = zeros(1,size(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.move_frames,:),1));
   trace_to_label(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.move_frames,:),2)==0) = 1;   
   pixellist = bwconncomp(trace_to_label);
%pixellist = bwconncomp(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.move_frames,:),2));
move_dropout_lengths{ll} = cellfun(@numel,pixellist.PixelIdxList);

   
fraction_missing_rest_nolong(ll) = (numel(find(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.rest_frames,:),2) == 0))-...
    sum(rest_dropout_lengths{ll}(rest_dropout_lengths{ll}>(long_thresh*300))))./...
       length_rest; 
   fraction_missing_move_nolong(ll) = (numel(find(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.move_frames,:),2) == 0))-...
    sum(move_dropout_lengths{ll}(move_dropout_lengths{ll}>(long_thresh*300))))./...
       length_move; 
   
end

set(0,'defaultfigurecolor',[1 1 1])

figure(555)
bar(fraction_missing_rest)
box off
set(gca,'XTick',1:numel(marker_names),'XTickLabels',marker_names)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Fraction of time missing')
title('Rest')
ylim([0 0.5])
xlim([0 numel(marker_names)+1])
print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerTracking_rest.png'))


figure(565)
bar(fraction_missing_rest_nolong)
box off
set(gca,'XTick',1:numel(marker_names),'XTickLabels',marker_names)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Fraction of time missing')
title('Rest')
ylim([0 0.5])
xlim([0 numel(marker_names)+1])



figure(556)
bar(fraction_missing_move)
box off
set(gca,'XTick',1:numel(marker_names),'XTickLabels',marker_names)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Fraction of time missing')
title('Move')

ylim([0 0.5])
xlim([0 numel(marker_names)+1])
print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerTracking_move.png'))



figure(566)
bar(fraction_missing_move_nolong)
box off
set(gca,'XTick',1:numel(marker_names),'XTickLabels',marker_names)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Fraction of time missing')
title('Move')

ylim([0 0.5])
xlim([0 numel(marker_names)+1])
print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerTracking_move_nolong.png'))



for ll = 1:numel(marker_names)

figure(557)

subplot_row = 4;
subplot_cols = ceil(numel(mocapstruct.markernames)./subplot_row);
 subplot_tight(subplot_row,subplot_cols,ll)
histvals = unique([0:50:600]);
    [n,x] = hist( rest_dropout_lengths{ll},histvals);
    bar(x,log10(n),'b');
%     hold on
%     [n2,x2] = hist(marker_velocity,0:0.1:5);
%     bar(x2,log10(n2),'r');
%     hold off
 title(mocapstruct.markernames{ll})
    xlim([0 600])
    xlabel('Marker Velociy [mm/frame]')
    ylabel('Fraction of bins')
  print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerDropoutDist_rest.png'))
  
  
figure(558)

subplot_row = 4;
subplot_cols = ceil(numel(mocapstruct.markernames)./subplot_row);
 subplot_tight(subplot_row,subplot_cols,ll)
histvals = unique([0:50:600,5000,10000]);
    [n,x] = hist( move_dropout_lengths{ll},histvals);
    bar(x,log10(n),'b');
%     hold on
%     [n2,x2] = hist(marker_velocity,0:0.1:5);
%     bar(x2,log10(n2),'r');
%     hold off
 title(mocapstruct.markernames{ll})
    xlim([0 600])
    xlabel('Marker Velociy [mm/frame]')
    ylabel('Fraction of bins')
  print('-dpng',strcat(mocapstruct.plotdirectory,'MarkerDropoutDist_move.png'))
end


end