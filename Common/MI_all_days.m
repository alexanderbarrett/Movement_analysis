
% 1. Tone on/of
% 2. Tone onset +/- 3 seconds
% 3. Tone offset +/-3 seconds
% 4. Extended tone trace (tone trace +/- 3 seconds)
%
% if shocking was done (only on day 3)
%
% 5. Shock on/off
% 6. Shock extended trace (Shock trace +/-2seconds)


% Run MI on all days


mouse = strcat('AM_imaging_3rd_set',filesep,'Mouse952');
num_days = 6;

for j =6:num_days;
    % if (j ~=3)
    strcat('day ',num2str(j))
    
    directory = strcat(mouse,filesep,'Session',num2str(j),filesep,'Analysis',filesep);
    
    
    traces = load(strcat(directory,'IC_traces_selected.mat'));
    raster = load(strcat(directory,'IC_events_selected.mat'));
    tone_file = importdata(strcat(directory,'Stimulus_traces.txt'),'\t');
        if(j~=3)
    freezing = load(strcat(directory,'Freezing_traces.mat'));
      else
               freezing = load(strcat(directory,'Freezing_traces.mat'));

    end
    
    
    
    
    framerate = 5;
    
    if(j~=3)
    freeze_mask = freezing.SpikeTraceData(5).Trace;
    else
            freeze_mask = freezing.SpikeTraceData(4).Trace;
    end
    
    tone_period = squeeze(tone_file(:,1));
    
    tone_onset = squeeze(tone_file(:,2));
    tone_offset = squeeze(tone_file(:,3));
    
    if (j==3)
        shock_correlates = squeeze(tone_file(:,6));
    end
    
    
    freezing_ind = find(freeze_mask == 1);
    moving_ind = find(freeze_mask == 0);
    
    numICs = min(size(raster.SpikeTraceData,2),size(traces.SpikeTraceData,2));
    length_trace = raster.SpikeTraceData(1).DataSize(1);
    
    
    raw_rasters_matrix = zeros(length_trace,numICs);
    raw_traces_matrix = zeros(traces.SpikeTraceData(1).DataSize(1),numICs);
    
    
    for i=1:size(raster.SpikeTraceData,2)
        temp = raster.SpikeTraceData(i).Trace;
        raw_rasters_matrix(1:length(temp),i) = temp;
    end
    
    
    num_shuff = 10000;
    
    
    if (j==3)
        [pvals_ret MI] = Mutual_Information_Shuffled(raw_rasters_matrix,shock_correlates,2,1,num_shuff);
        [pvals_ret MI] = Mutual_Information_Shuffled(raw_rasters_matrix,[freeze_mask tone_onset ...
            tone_offset tone_period shock_correlates],2,1,num_shuff);
        
        pvals_freeze = pvals_ret(:,1); MI_freeze = MI(:,1);
        save(strcat(directory,'mi_test_freeze_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
        pvals_freeze = pvals_ret(:,2); MI_freeze = MI(:,2);
        save(strcat(directory,'mi_test_toneon_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
        pvals_freeze = pvals_ret(:,3); MI_freeze = MI(:,3);
        save(strcat(directory,'mi_test_toneoff_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
        pvals_freeze = pvals_ret(:,4); MI_freeze = MI(:,4);
        save(strcat(directory,'mi_test_toneperiod_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
        pvals_freeze = pvals_ret(:,5); MI_freeze = MI(:,5);
        save(strcat(directory,'mi_test_shock_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
    else
        
        [pvals_ret MI] = Mutual_Information_Shuffled(raw_rasters_matrix,[freeze_mask tone_onset...
            tone_offset tone_period],2,1,num_shuff);
        
        pvals_freeze = pvals_ret(:,1); MI_freeze = MI(:,1);
        save(strcat(directory,'mi_test_freeze_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
        pvals_freeze = pvals_ret(:,2); MI_freeze = MI(:,2);
        save(strcat(directory,'mi_test_toneon_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
        pvals_freeze = pvals_ret(:,3); MI_freeze = MI(:,3);
        save(strcat(directory,'mi_test_toneoff_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
        pvals_freeze = pvals_ret(:,4); MI_freeze = MI(:,4);
        save(strcat(directory,'mi_test_toneperiod_',num2str(num_shuff),'.mat'),'pvals_freeze','MI_freeze')
        
    end
    
    
    
    
end