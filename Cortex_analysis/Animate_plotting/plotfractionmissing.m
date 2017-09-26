function plotfractionmissing(mocapstruct)
marker_names = fieldnames(mocapstruct.markers_preproc);

fraction_missing_rest = zeros(1,numel(marker_names));
fraction_missing_move = zeros(1,numel(marker_names));

length_rest = numel(mocapstruct.rest_frames);
length_move = numel(mocapstruct.move_frames);
%size_mocap = size( mocapstruct.markers_preproc.(marker_names{1}),1);
for ll = 1:numel(marker_names)
   fraction_missing_rest(ll) = numel(find(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.rest_frames,:),2) == 0))./...
       length_rest; 
   fraction_missing_move(ll) = numel(find(sum(mocapstruct.markers_preproc.(marker_names{ll})(mocapstruct.move_frames,:),2) == 0))./...
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
end