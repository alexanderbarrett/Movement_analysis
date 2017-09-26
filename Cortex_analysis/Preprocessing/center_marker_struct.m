function [markers_centered] = center_marker_struct(markers,fps)


markers_centered = markers;
markernames = fieldnames(markers);

%% remove high frequency noise (>10 Hz) in the spinal trace

dH = designfilt('lowpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 10/(fps/2), ...
    'DesignMethod', 'butter');
[f1,f2] = tf(dH);

marker_M_lowpass =zeros(size(markers.SpineM));
markers.SpineM(isnan(markers.SpineM)) = 0;
for mk = 1:3        
marker_M_lowpass(:,mk) = filtfilt(f1,f2,...
            markers.SpineM(:,mk));
end

%% subtract from other markers
for jj = 1:numel(markernames)
    markers_centered.(markernames{jj}) = markers_centered.(markernames{jj})-marker_M_lowpass;
end

end