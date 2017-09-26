function plot_marker_lever_xcorr(markers,lever_thresholded,analog)


marker_names = fieldnames(markers);
num_markers = numel(marker_names);
analog_names = fieldnames(analog);

trace_length = size(markers.(marker_names{1}),1);

index_xcorr = 1:numel(lever_thresholded);

figure(33333)
[vals,lags] = xcorr(lever_thresholded(index_xcorr),markers.ElbowR(index_xcorr,3),2000,'coeff');
plot(lags,vals,'b')
[~,maxind] = max(vals);
lags(maxind)
hold on
[vals,lags] = xcorr(lever_thresholded(index_xcorr),markers.HeadF(index_xcorr,3),2000,'coeff');
plot(lags,vals,'g')
[~,maxind] = max(vals);
lags(maxind)
[vals,lags] = xcorr(lever_thresholded(index_xcorr),markers.ElbowL(index_xcorr,3),2000,'coeff');
plot(lags,vals,'r')
[~,maxind] = max(vals);
lags(maxind)
hold off
xlabel('Time lag between analog and R elbow')
ylabel('Cross correlation Coefficient')
box off

figure(2344)
plot(lever_thresholded,'b')
hold on
plot(markers.ElbowR(:,3)*0.02,'r');
[~,maxind] = max(vals);
plot(circshift(lever_thresholded,[2 -lags(maxind)]),'g')
hold off
