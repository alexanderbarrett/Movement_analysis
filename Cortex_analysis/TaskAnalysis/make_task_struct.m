function taskstruct = make_task_struct(analog,lever_thresholded,fps)
 

%% start from scratch and get thresholded and bout labeled licks, speakers and taps

    %% get lick onsets
    speaker_on = find( abs(analog.SPEAKER)>0.3);
    speaker_on_times = zeros(1,numel(analog.SPEAKER));
    speaker_on_times(speaker_on) = 1;
    %length of one cycle
    speaker_on_times = conv(speaker_on_times,ones(1,5)./5,'same');
    speaker_on_times( speaker_on_times >0) = 1;    
    [speaker_labels,num_tones] = bwlabel(speaker_on_times);
    
    %% get lick onsets and put into lick bouts
    lick_on = find( abs(analog.LICK_SENSOR)>0.3);
    lick_on_times = zeros(1,numel(analog.LICK_SENSOR));
    lick_on_times(lick_on) = 1;
    lick_on_times = conv(lick_on_times,ones(1,fps/3)./(fps/3));
    lick_on_times( lick_on_times >0) = 1;
    
    [lick_labels,num_licks] = bwlabel(lick_on_times);
    
     %% resample the lever and plot
    resampled_lever_thresh = [0 diff(lever_thresholded)];
    resampled_lever_thresh(resampled_lever_thresh>0.05) = 1;
    resampled_lever_thresh(resampled_lever_thresh<-0.05) = -1;
    resampled_lever_thresh(intersect(find(resampled_lever_thresh<0.1),...
        find(resampled_lever_thresh>-0.1))) = 0;
    
    tap_times =find(resampled_lever_thresh==1);
    
    %% get lick and speaker onsets for bouts
     speaker_onsets = zeros(1,num_tones);
    lick_onsets = zeros(1,num_licks);
    
    %% get the first lick and tonie within a bout
    for ll = 1:num_tones
        speaker_onsets(ll) = floor(find(speaker_labels == ll,1,'first'));
    end
    
    for ll = 1:num_licks
        lick_onsets(ll) = floor(find(lick_labels == ll,1,'first'));
    end
    
    
    
    %% plot the data
    figure(1234)
    plot(analog.LICK_SENSOR)
    hold on
    plot(lick_on_times,'r')
    plot(speaker_on_times,'g')
    plot(lever_thresholded,'c')
    
       figure(324)
%    plot(x_axis,resampled_lever_thresh,'k')
    
    hold on
    %plot(x_axis_analog,lever_thresholded,'g')
    hold off
    
    
    %% get time around
    time_pre_lick = 0;
    time_post_lick = 3;
    timing_lick = -time_pre_lick*fps:time_post_lick*fps;
    
    
    
    %% get the events within the speaker times
    good_lick_times = reshape(bsxfun(@plus,speaker_onsets',timing_lick),1,[]);
    
    %% get the incorrect lick onsets and the correct
    [success_licks,ind1,ind2] = intersect(lick_onsets,good_lick_times);
    lick_onsets_bad = lick_onsets;
    lick_onsets_bad(ind1) = [];

good_licks = intersect(lick_on_times,good_lick_times );
bad_licks = setxor(lick_on_times,good_licks );



%% now sort good and bad trials
intertrial_frames = ceil(fps*1.2);
    trial_diff = diff(tap_times);
    num_trials = 1;
    trial_tap_inds = cell(1,1); 
    trial_tap_inds{1} = 1;
    
    for mm = 1:numel(tap_times)-1;
    if (trial_diff(mm)<intertrial_frames)
        trial_tap_inds{num_trials} = cat(1,trial_tap_inds{num_trials},mm+1);
    else
       num_trials=num_trials+1;
       trial_tap_inds{num_trials} = mm+1;
    end    
    end
    
    
    %% get tap times
    tap1_times = zeros(1,num_trials);
    tap2_times = zeros(1,num_trials);
    tapn_times = cell(1,num_trials);
    good_trials = zeros(num_trials);
    for mm = 1:num_trials
        tap1_times(mm) = tap_times(trial_tap_inds{mm}(1));
        if numel(tap_times(trial_tap_inds{mm}))>1
            tap2_times(mm)  = tap_times(trial_tap_inds{mm}(2));
        else
          tap2_times(mm) = nan;
        end
            if numel(tap_times(trial_tap_inds{mm}))>2
    tapn_times{mm}  = tap_times(trial_tap_inds{mm}(2:end));
            end
    end
    
    %% get good trials
    speaker_levers = zeros(1,numel(speaker_onsets));
     for kk = 1:numel(speaker_onsets)
        difference_here = tap_times-speaker_onsets(kk);
        [diffval,ind_min] = find(difference_here<0,1,'last'); %find the tap preceeding
       difference_here(ind_min);
      
        %if there is at most 0.5 s between tap and tone
        if (abs(difference_here(ind_min))<0.25*fps)
        speaker_levers(kk) = tap_times(ind_min);
        else
             speaker_levers(kk) = nan;
        end
     end
    
     [~,good_trials,~] = intersect(tap2_times,speaker_levers);
     
     taskstruct.tap1_times = tap1_times;
    taskstruct.tap2_times = tap2_times;
    taskstruct.tapn_times = tapn_times;
    taskstruct.good_trials = good_trials;
    taskstruct.num_trials = num_trials;
        taskstruct.bad_trials = setxor(1:num_trials,good_trials);
        taskstruct.good_lick_times = good_licks;
        taskstruct.bad_lick_times = bad_licks;
        taskstruct.resampled_lever_thresh = resampled_lever_thresh;

end
    