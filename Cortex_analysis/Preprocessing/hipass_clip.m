function [marker_clipped,clipped_index] = hipass_clip(markers,bad_frames,params)


dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 0.3/(params.fps/2), ...
    'DesignMethod', 'butter');
[f1_hipass,f2_hipass] = tf(dHipass);



marker_clipped = struct();
markernames = fieldnames(markers);
marker_clipped = markers;
%cluster_size_dist = cell(1,number_of_clust(mmm));

inst_label = ones(1,numel(markers.(markernames{1})(:,3)));
inst_label(bad_frames) = 0;
clipped_index = 1:numel(inst_label);
pixellist = bwconncomp(inst_label);
for mk = 1:numel(markernames)
    fprintf('Clipping marker %s \n',markernames{mk})
    for lk =1:pixellist.NumObjects
        if (numel(pixellist.PixelIdxList{lk})>30)
        for sk = 1:3
            marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk) = ...
                filtfilt(f1_hipass,f2_hipass,marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk)-...
                mean(marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk)));
        end
        else
            inst_label(pixellist.PixelIdxList{lk}) = 0;
        end
        
    end
    
    
    %% also remove the bad frames
end
for mk = 1:numel(markernames)
      marker_clipped.(markernames{mk})(inst_label == 0,:) = [];
end
    clipped_index(inst_label==0) = [];


end