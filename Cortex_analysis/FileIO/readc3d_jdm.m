
function [markers,analog,resample_analog,lever_thresholded] = readc3d_jdm(filepath,fps,analog_fps)

acq = btkReadAcquisition(filepath);

markers = btkGetMarkers(acq);
analog = btkGetAnalogs(acq);

btkCloseAcquisition(acq);



%% get frame rate
fid = fopen(filepath, 'r');
fseek(fid,0,'cof');
frame_rate=fread(fid,1,'float32');
fclose(fid);


%% Process the analog data

w0 = 60/(analog_fps/2); %filter fundamental frequency
bw = w0/50; %Q=35
[b,a] = iirnotch(w0,bw);
filtered_lever = filtfilt(b,a,analog.LEVER);

lever_threshold = 0.6;

lever_thresholded = zeros(1,numel(filtered_lever));
lever_thresholded(filtered_lever>lever_threshold) = 1;
lever_thresholded(filtered_lever<-lever_threshold) = -1;

figure(111)
plot(filtered_lever)
hold on
plot(0.6*ones(1,numel(filtered_lever)),'r')
plot(analog.LICK_SENSOR,'g')
plot(lever_thresholded,'k')
plot(analog.SPEAKER,'y')
hold off


resample_analog = resample([lever_thresholded' analog.SPEAKER analog.WATER_VALVE analog.LICK_SENSOR analog.GROUND],1,20);
resample_analog = bsxfun(@minus,resample_analog,mean(resample_analog,1));
end