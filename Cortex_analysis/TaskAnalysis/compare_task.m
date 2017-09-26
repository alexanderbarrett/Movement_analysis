function compare_task(mocapstruct,mocapstruct_2)
fps = params_task.fps
analog_fps = fps;

marker_frame_length = numel(lever_thresholded);
x_axis = 0:1/fps:(marker_frame_length-1)./fps;
x_axis_analog = 0:1/analog_fps:(numel(lever_thresholded)-1)./analog_fps;


    taskstruct = make_task_struct(analog,lever_thresholded,fps);

%     

  
    %% plot lever and speaker
        figure(1239)
  %  plot(x_axis,0.001*xcorr_val(lagval>-1))
    hold on
    plot(x_axis_analog,analog.SPEAKER,'y')
    plot(x_axis_analog,lever_thresholded,'g')
    % plot(x_axis_analog,speaker_on_times,'r')
          plot(x_axis_analog,taskstruct.resampled_lever_thresh,'b')

    hold off
   
    

    
    %% get the lever onsets of unsuccessful periods
%     lever_onsets_bad = zeros(1,numel(lick_onsets_bad));
%     for kk = 1:numel(lick_onsets_bad)
%         difference_here = find(resampled_lever_thresh==1)-lick_onsets_bad(kk);
%         tap_times =find(resampled_lever_thresh==1);
%         [~,ind_min] = find(difference_here<0,1,'last'); %find the tap preceeding
%         if (numel(ind_min))
%         lever_onsets_bad(kk) = tap_times(ind_min);
%         end
%     end
%     
%     %% get the lever onsets of successful periods
%     lever_onsets_good = zeros(1,numel(speaker_onsets));
%             tap_times =find(resampled_lever_thresh==1);
% 
%     for kk = 1:numel(speaker_onsets)
%         difference_here = find(resampled_lever_thresh==1)-speaker_onsets(kk);
%         [diffval,ind_min] = find(difference_here<0,1,'last'); %find the tap preceeding
%        difference_here(ind_min)
%       
%         %if there is at most 0.5 s between tap and tone
%         if (abs(difference_here(ind_min))<0.25*fps)
%         lever_onsets_good(kk) = tap_times(ind_min);
%         else
%              lever_onsets_good(kk) = nan;
%         end
%     end
%     lever_onsets_good(isnan(lever_onsets_good)) = [];
  
    %marker_velocity(13,lever_template,4);
    
    
lever_onsets_good = taskstruct.tap2_times(taskstruct.good_trials);
    lever_onsets_bad= taskstruct.tap1_times(taskstruct.bad_trials);

    
    
    %% get time pre post for plotting
    time_pre = 2;
    time_post = 2;
    timing = floor(-time_pre*fps:time_post*fps);
    %%
    time_pre_template = 1.3;
    time_post_template = 0.2;
    timing_template = floor(-time_pre_template*fps:time_post_template*fps);
    
    %% plot relative to the speaker onset
    reach_times = bsxfun(@plus,lever_onsets_good',timing);
        success_lever_times = floor(bsxfun(@plus,lever_onsets_good',timing));

    
    unsuccess_reach_times = floor(bsxfun(@plus,lever_onsets_bad',timing));
    unsuccess_lever_times = floor(bsxfun(@plus,lever_onsets_bad',timing));
    
    unsuccess_reach_times(unsuccess_reach_times<1) = 1;
    unsuccess_lever_times(unsuccess_lever_times<1) = 1;
    
  %  unsuccess_lever_times(unsuccess_lever_times<1) = 1;
  %  unsuccess_lever_times(unsuccess_lever_times<1) = 1;
    
    
    %unsuccess_reach_times(unsuccess_reach_times<1) = 1;
    markernames = fieldnames(markers_preproc);
    for ll = 1:numel(markernames)
        markers_preproc.(markernames{ll})(find(sum(markers_preproc.(markernames{ll}),2)==0),:) = nan;
    end
    
    
    lhand = markers_preproc.ElbowL(:,3);
    lhandx = markers_preproc.ElbowL(:,1);
    lhandy = markers_preproc.ElbowL(:,2);
    lhand_vel = conv(marker_velocity(13,:,3),ones(1,20)./20,'same');
    
    rhand = markers_preproc.ElbowR(:,3);
    fhead = markers_preproc.HeadF(:,3);
    rspine = markers_preproc.SpineF(:,3);
    
    
    lever_template = bsxfun(@plus,lever_onsets_good',timing_template);
    lhandz_template = (lhand(lever_template(1,:))-mean(lhand(lever_template(1,:))))./std(lhand(lever_template(1,:)));
    [xcorr_val,lagval] = xcorr(lhand,lhandz_template);
    

    
    
    %% visualize the successful taps
            trials_plot = intersect(1:1:50,1:size(success_lever_times));

    fighere = figure(233)
    set(fighere,'color','w')
    subplot(1,3,1)
    plot(lhand(success_lever_times(trials_plot,:))')
       hold on
    plot(5*sum(taskstruct.resampled_lever_thresh(success_lever_times(trials_plot,:))',2))
    hold off
    
    box off
    xlabel('Time (frames)')
    ylabel('z position (mm)')
    title('left hand')
    
    subplot(1,3,2)
    % plot(lhandx(success_lever_times(10:12,:))')
    plot(fhead(success_lever_times(trials_plot,:))')
    hold on
    plot(5*sum(taskstruct.resampled_lever_thresh(success_lever_times(trials_plot,:))',2))
    hold off
    
    xlabel('Time (frames)')
    ylabel('z position (mm)')
    title('front of head')
    box off
    
    subplot(1,3,3)
    plot(rhand(success_lever_times(trials_plot,:))')
       hold on
    plot(5*sum(taskstruct.resampled_lever_thresh(success_lever_times(trials_plot,:))',2))
    hold off
    
    xlabel('Time (frames)')
    ylabel('z position (mm)')
    title('right hand')
    box off
    if (plotdirectory ~= 0)
    print('-depsc',strcat(plotdirectory,'success_levertaps.eps'))
    print('-dpng',strcat(plotdirectory,'success_levertaps.png'))
    end
    
        trials_plot = intersect(1:5:250,1:size(reach_times));

    figure(234)
    subplot(1,3,1)
    plot(lhand(reach_times(trials_plot,:))')
    
    subplot(1,3,2)
    plot(fhead(reach_times(trials_plot,:))')
    
    subplot(1,3,3)
    plot(rspine(reach_times(trials_plot,:))')
    
    %             subplot(1,3,3)
    %           plot(fhead(reach_times)')
    %
    % hold on
    %plot(speaker_on_times,'r')
    
    trials_plot = intersect(1:1:50,1:size(unsuccess_lever_times));
    
    fighere = figure(235)
    set(fighere,'color','w')
    subplot(1,3,1)
    plot(lhand(unsuccess_lever_times(trials_plot,:))')
    box off
    xlabel('Time (frames)')
    ylabel('z position (mm)')
    title('left hand')
    
    subplot(1,3,2)
    % plot(lhandx(success_lever_times(10:12,:))')
    plot(fhead(unsuccess_lever_times(trials_plot,:))')
    xlabel('Time (frames)')
    ylabel('z position (mm)')
    title('front of head')
    box off
    
    subplot(1,3,3)
    plot(rhand(unsuccess_lever_times(trials_plot,:))')
    xlabel('Time (frames)')
    ylabel('z position (mm)')
    title('right hand')
    box off
    if (plotdirectory ~= 0)
    print('-depsc',strcat(plotdirectory,'unsuccess_levertaps.eps'))
    print('-dpng',strcat(plotdirectory,'unsuccess_levertaps.png'))
    end
    %% make a task visualization
    
    frame_inds = reshape(success_lever_times(trials_plot,:)',[],1);
    frame_inds = frame_inds(1:10:500);
    h=figure(444);
    view([90, 12]);
    M=[];
    
    if (plotdirectory ~= 0)
      v = VideoWriter(strcat(plotdirectory,'success_movie.mp4'),'MPEG-4');
     open(v)
    end
        [markers_aligned_here]=  center_marker_struct(markers_preproc,fps);
       %  [~,markers_aligned_here]= align_hands_elbows_task(markers_preproc,fps);
        % center_marker_struct(markers_preproc,fps)
  animate_markers_clusteraligned(markers_aligned_here,markers_preproc,...
      frame_inds,marker_names,markercolor,links,M,1,'1',h,1)
      if (plotdirectory ~= 0)
  close(v)
    end
