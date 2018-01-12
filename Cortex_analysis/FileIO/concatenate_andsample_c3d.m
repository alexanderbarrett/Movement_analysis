
function [marker_agg,analog_agg,resample_analog_agg,lever_thresholded_agg,filestartpts] = concatenate_c3d(filepath_array,fps,analog_fps,desired_length,chunksize)
marker_agg = [];
analog_agg = [];
resample_analog_agg = [];
lever_thresholded_agg = [];
filestartpts = zeros(2,numel(filepath_array));

correct_analog_time = 1;

for ll = 1:numel(filepath_array)
acq = btkReadAcquisition(filepath_array{ll});

markers = btkGetMarkers(acq);
analog = btkGetAnalogs(acq);

btkCloseAcquisition(acq);

fnames = fieldnames(analog);

%% get frame rate
fid = fopen(filepath_array{ll}, 'r');
fseek(fid,0,'cof');
frame_rate=fread(fid,1,'float32');
fclose(fid);


%% Process the analog data

w0 = 60/(analog_fps/2); %filter fundamental frequency
bw = w0/50; %Q=35
[b,a] = iirnotch(w0,bw);

   fn_markers = fieldnames(markers);    
marker_number = (numel(markers.(fn_markers{1})(:,1)));

  fn_analog = fieldnames(analog);   

if (numel(fn_analog) == 0)
    analog = struct('GROUND',zeros((marker_number),1),...
        'LEVER',zeros((marker_number),1),...
        'WATERVALVE',zeros((marker_number),1),...
        'SPEAKER',zeros((marker_number),1),...
        'LICK_SENSOR',zeros((marker_number),1));
    

%analog.LICK_SENSOR = zeros(1,numel(marker_number));
lever_thresholded = zeros(1,(marker_number));

else
    


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


%% supplement the length if there is a mismatch



%% there is an offset between the analog and the true, correct for it here
xcorr_chunk = 10^5;
lever_shifted = lever_thresholded;
analog_shifted = analog;
analog_number = numel(analog.(fn_analog{1}));

if analog_number <= marker_number
for jk =1:xcorr_chunk:numel(lever_thresholded)
    
    
    ind_start = jk;
    ind_finish = jk+xcorr_chunk;
    if (ind_finish>numel(lever_thresholded))
        ind_finish = numel(lever_thresholded);
    end
    index_xcorr = ind_start:ind_finish;
    
    if numel(find(lever_thresholded(index_xcorr)))>100
figure(333)
[vals,lags] = xcorr(lever_thresholded(index_xcorr),markers.ElbowR(index_xcorr,3),2000,'coeff');
plot(lags,vals,'b')
[~,maxind_1] = max(vals);

hold on
[vals,lags] = xcorr(lever_thresholded(index_xcorr),markers.HeadF(index_xcorr,3),2000,'coeff');
plot(lags,vals,'g')
[~,maxind_2] = max(vals);

[vals,lags] = xcorr(lever_thresholded(index_xcorr),markers.ElbowL(index_xcorr,3),2000,'coeff');
plot(lags,vals,'r')
[~,maxind_3] = max(vals);



lag_shift = floor(mean([lags(maxind_3) lags(maxind_2) lags(maxind_1)]));
hold off

    else
        lag_shift = 400;
    end

lever_shifted(index_xcorr) = circshift(lever_shifted(index_xcorr),-lag_shift,2);
num_xcorr = numel(index_xcorr);
if (num_xcorr-lag_shift >0)
lever_shifted(index_xcorr(num_xcorr-lag_shift:num_xcorr)) = 0;
end
    
    
for mm = 1:numel(fn_analog)
    analog_shifted.(fn_analog{mm})(index_xcorr) = circshift(analog_shifted.(fn_analog{mm})(index_xcorr), -lag_shift);
   analog_shifted.(fn_analog{mm})(index_xcorr(num_xcorr-lag_shift:num_xcorr)) = 0;

end


end
end
analog = analog_shifted;
lever_thresholded = lever_shifted;

end
  fn_analog = fieldnames(analog);   

analog_number = numel(analog.(fn_analog{1}));


resample_analog = resample([lever_thresholded' analog.SPEAKER analog.(fn_analog{3}) analog.LICK_SENSOR analog.GROUND],1,1);
resample_analog = bsxfun(@minus,resample_analog,mean(resample_analog,1));


if analog_number < marker_number
    fprintf('Mismatch for marker and analaog lengths for file number %f %s \n',ll,filepath_array{ll});
   for jj = 1:numel(fn_analog)
       analog.(fn_analog{jj}) = cat(1,analog.(fn_analog{jj}),zeros(marker_number-analog_number,1));
   end
   resample_analog = cat(1,resample_analog,zeros(marker_number-analog_number,size(resample_analog,2)));
     lever_thresholded = cat(2,lever_thresholded,zeros(size(lever_thresholded,2),marker_number-analog_number));

end





if (ll == 1)
      fn_markers = fieldnames(markers);    
      
      virtualmarkers = [];
      realmarkers = markers;
      for mm = 1:numel(fn_markers)
          if (strcmp(fn_markers{mm}(1),'V') || strcmp(fn_markers{mm},'Head_between'))
              virtualmarkers = cat(1,virtualmarkers,mm);
              realmarkers = rmfield(realmarkers,fn_markers{mm});
          end
      end
          good_markers = setxor(1:numel(fn_markers),virtualmarkers);
      
    marker_agg = realmarkers;
analog_agg = analog;
resample_analog_agg = resample_analog;
lever_thresholded_agg = lever_thresholded;

filestartpts(1,ll) = 1;
filestartpts(2,ll) = size(markers.(fn_markers{1}),1);

else
    %% concatenate analog
    fn_analog = fieldnames(analog);    
    %analog_agg=struct();
   for jj = 1:numel(fn_analog)
       analog_agg.(fn_analog{jj}) = cat(1,analog_agg.(fn_analog{jj}),analog.(fn_analog{jj}));
   end
   
      %% concatenate analog resampled
      resample_analog_agg = cat(1,resample_analog_agg,resample_analog);
%    fn_analog = fieldnames(resample_analog);    
%    for jj = 1:numel(fn_analog)
%        resample_analog_agg.(fn_analog{jj}) = cat(1,resample_analog_agg.(fn_analog{jj}),resample_analog.(fn_analog{jj}));
%    end
   
      %% concatenate markers
      filestartpts(1,ll) = size(marker_agg.(fn_markers{1}),1)+1;

   fn_markers = fieldnames(marker_agg);    
   for jj = 1:numel(fn_markers)
       marker_agg.(fn_markers{jj}) = cat(1,marker_agg.(fn_markers{jj}),markers.(fn_markers{jj}));
   end
   
   %% concatenate thresholded lever
   lever_thresholded_agg = cat(2,lever_thresholded_agg,lever_thresholded);
    
   
filestartpts(2,ll) = size(marker_agg.(fn_markers{1}),1);
   
end
end


markers = marker_agg;
analog = analog_agg;
resample_analog = resample_analog_agg;
lever_thresholded = lever_thresholded_agg; 

   marker_names = fieldnames(markers);
    num_markers = numel(marker_names);
    analog_names = fieldnames(analog);
    
    
    trace_length = size(markers.(marker_names{1}),1);
    
    %% reinitialize
        for lk = 1:numel(marker_names)
            marker_agg.(marker_names{lk}) =  [];
        end
        for lk=1:numel(analog_names)
            analog_agg.(analog_names{lk}) =  [];
        end
  
     lever_thresholded_agg=[];
     resample_analog_agg = [];
%% randomly take chunks from data to examine a subset
if (desired_length<size(lever_thresholded))
    numchunks = floor(desired_length./chunksize);
    gaplength = trace_length-desired_length;
    gapsize = floor(gaplength./numchunks);
    chunkedstarts = 1:(chunksize+gapsize):(numchunks-1)*(chunksize+gapsize);
    frames_to_use = sort(reshape(bsxfun(@plus,1:chunksize,chunkedstarts'),1,[]),'Ascend');
    frames_to_use=unique(frames_to_use);
    
    lever_thresholded = lever_thresholded(frames_to_use);
    resample_analog = resample_analog(frames_to_use,:) ;
    
     %% this is redundant but I'm too lazy to change it
    for lk = 1:numel(marker_names)
        markers.(marker_names{lk}) =  markers.(marker_names{lk})(frames_to_use,:);
    end
    for lk=1:numel(analog_names)
        analog.(analog_names{lk}) =  analog.(analog_names{lk})(frames_to_use);
    end
end
    

    
   
    
   
    %% aggregate
    lever_thresholded_agg = cat(1,lever_thresholded_agg,lever_thresholded);
    resample_analog_agg = cat(1,resample_analog_agg,resample_analog);
    
    for lk = 1:numel(marker_names)
        marker_agg.(marker_names{lk}) =  cat(1,marker_agg.(marker_names{lk}),markers.(marker_names{lk}));
    end
    for lk=1:numel(analog_names)
        analog_agg.(analog_names{lk}) =  cat(1,analog_agg.(analog_names{lk}),analog.(analog_names{lk}));
    end
    
end







