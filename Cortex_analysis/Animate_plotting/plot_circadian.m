function plot_circadian(mocapstruct)

figure(5660)
plot(mocapstruct.markers_preproc.HeadF(:,3))
set(gca,'XTick',mocapstruct.filestartpts(1,:),'XTickLabel',mocapstruct.mocapfiletimes)