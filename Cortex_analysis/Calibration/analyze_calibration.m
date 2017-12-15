stabilityfile = {'X:\Jesse\Data\Motionanalysis_captures\Calibration_recordings\20171027\calibration_recording_stability1.c3d'};


desired_length = 8*10^6;
chunksize = 300*fps;
[markers,analog,resample_analog,lever_thresholded,filestartpts] = concatenate_andsample_c3d(stabilityfile,300,300,...
    desired_length,chunksize);



markerfields = fieldnames(markers); 
numinds = 3;
std_vals = zeros(numel(markerfields),numinds);
for jj = 1:numel(markerfields)
    for mm = 1:3
        std_vals(jj,mm) = std(markers.(markerfields{jj})(:,mm));
   % [x,n] = hist(bsxfun(@minus,markers.(markerfields{jj})(:,mm),nanmean(markers.(markerfields{jj})(:,mm),2)),-1:0.001:1);
   % plot(n,x);
    end
end
    
    mean(std_vals,1)
figure(44)
plot(0:1./300:(size(markers.M1,1)-1)./300,markers.M1(:,1))
xlabel('Time (s)')
ylabel('X position (mm)')
title('X position of a static marker')
box off