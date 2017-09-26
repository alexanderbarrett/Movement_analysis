

[markers,markers_aligned] = align_hands_elbows(markers,fps);

marker_frame_length = size(markers.(marker_names{1}),1);

%%initialize
markers_preproc = markers;
markers_preproc_aligned = markers_aligned;



%% other interpolation etc.
for ll = 1:numel(marker_names)
    fprintf('starting interpolation for markers, over 30 frames (100 ms) maximum %f \n',ll)
    % for jkj = 1:3 %need to do over all x y z simul and add spike correction
    [temp,fake_frames] = markershortinterp((markers_preproc.(marker_names{ll})),30,5);
    fakedata.(marker_names{ll}){1} = fake_frames;
    % if numel(temp)
    markers_preproc.(marker_names{ll}) = temp;
    %  end
    clear fake_frames
    
    fprintf('interpolating aligned markers \n')
    [temp,~] =markershortinterp((markers_preproc_aligned.(marker_names{ll})),100,5);
    % if numel(temp)
    markers_preproc_aligned.(marker_names{ll}) = temp;
    %  end
end








%'big-data' features
% get relative marker positions to one another (x,y,z)
%delta_markers = zeros(num_markers,num_markers,marker_frame_length,4);
marker_velocity = zeros(num_markers,marker_frame_length,4);
marker_position = zeros(num_markers,marker_frame_length,3);
abs_velocity_antialiased = zeros(num_markers,marker_frame_length);
%marker_position_hipass = zeros(num_markers,marker_frame_length,3);


dH = designfilt('lowpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 60/(fps/2), ...
    'DesignMethod', 'butter');
[f1,f2] = tf(dH);

%delta_markers_reshaped = [];
fprintf('getting velocities \n')
for ll = 1:numel(marker_names)
    marker_position(ll,:,1:3) = markers_preproc.(marker_names{ll});
    for mk = 1:3
%         marker_position_hipass(ll,:,mk) = filtfilt(f1_hipass,f2_hipass,...
%             markers_preproc.(marker_names{ll})(:,mk));
        %
        %         tester = filtfilt(f1_hipass,f2_hipass,...
        %             markers_preproc.(marker_names{ll})(1:500,mk));
    end
    
    marker_velocity(ll,2:(end),1:3) = diff(markers_preproc.(marker_names{ll}));
    marker_velocity(ll,1,1:3) = marker_velocity(ll,2,1:3);
    marker_velocity(ll,:,4) = sqrt(sum((squeeze( marker_velocity(ll,:,1:3))).^2,2));
    abs_velocity_antialiased(ll,:) =  filtfilt(f1,f2, marker_velocity(ll,:,4));
    
    
    for lk = (ll+1):num_markers
        %delta_markers(ll,lk,:,1:3) =   markers_preproc_aligned.(marker_names{ll})-markers_preproc.(marker_names{lk});
        distance_here =   (markers_preproc.(marker_names{ll})-markers_preproc.(marker_names{lk}));
        %  delta_markers(ll,lk,:,4) =sqrt(sum((distance_here).^2,2));
        
        
        % delta_markers_reshaped = cat(2,  delta_markers_reshaped,markers_preproc.(marker_names{ll})-markers_preproc.(marker_names{lk}));
        
    end
end




%get aggregate feature matrix
%% simple bad frame detector
fprintf('finding bad frames \n')
tic


bad_frame_surround = fps./6;
bad_frames_agg = cell(1,num_markers);
bad_frame_velthresh = 20;
for mm = 1:num_markers
    bad_1 = find(squeeze(marker_velocity(mm, :,4))>bad_frame_velthresh);
    [i3] = ind2sub(size(squeeze(marker_velocity(mm, :,4))),bad_1);
    clear bad_1
    
    % find zeros - no marker found
    bad_2 = find(squeeze(sum(marker_position(mm,:, :),3))==0);
    [i4] = ind2sub(size(squeeze(sum(marker_position(mm,:, :),3))),bad_2);
    clear bad_2
    
    if numel(i4)
        bad_frames  = (cat(2,(i3),(i4)));
    else
        bad_frames = i3;
    end
    clear i3 i4
    
    bad_frames = unique((bad_frames));
    numframes = size(marker_velocity,2);
    gap_base = zeros(1,numframes );
    gap_base(bad_frames) = 1;
    
    gap_base_conv = conv(gap_base,ones(1,bad_frame_surround)./bad_frame_surround,'same');
    %gap_base_conv(gap_base_conv>0) = 1;
    bad_frames = find(gap_base_conv >0);
    
    bad_frames = reshape(bad_frames,1,[]);
    bad_frames = unique(bad_frames);
    bad_frames(bad_frames<1)=1;
    bad_frames(bad_frames>numframes ) = numframes ;
    
    bad_frames_agg{mm} = bad_frames;
end
clear marker_position

toc




cluster_here = [2];

    %head, L arm R arm, head and hips, hips, full
    cluster_markersets = {[1:3],[1:10,17,18],[4,7,11:13],[4,7,14:16],[6,8,9,18,19],[6,10,17,20],[6,8,9,10,17:20],[1:numel(marker_names)]};
    % cluster_markersets = {[1,2,3],[6,11,12,13],[5,9,10,17,18],[6,14,15,16],[4,5,6,7,8],[1:numel(marker_names)]};
    subcluster_names = {'head','axial','Larm','Rarm','L leg','R leg','Both legs','Global'};
    
    % subcluster_names = {'head','Larm','Hips and knees','Rarm','Spine','Global'};
    %        cluster_markersets = {[1:numel(marker_names)],[1,2,3],[3,4,5,6,7],[5,6,7,8],[1,2,3,9,10]};
    number_of_clust = [20,50,100,20,100,100,100,200];


number_of_subclusters = numel(cluster_markersets);

cluster_objects = cell(1,number_of_subclusters);

clustering_inds_agg = cell(1,number_of_subclusters);
badframe_inds_agg = cell(1,number_of_subclusters);

for mmm=[1:number_of_subclusters]
    temp = zeros(1,trace_length );
    cluster_marker_inds = cluster_markersets{mmm};
    for ll = cluster_marker_inds
        temp(bad_frames_agg{ll}) = 1;
        %clustering_inds{mmm} = cat(2,clustering_inds{mmm},bad_frames_agg{ll});
    end
    %for aligned pose, also get rid of the
    % temp(bad_frames_agg{5}) = 1;
    
    clustering_inds_agg{mmm} = find(temp == 0);
    badframe_inds_agg{mmm} = find(temp == 1);
end







%% plot marker characteristics to get a feel for overall velocities, range of motion, etc. 
params.fps= fps;
%find_rest_frames = (markers_preproc,marker_velocity,
plot_marker_characteristics(markers_preproc,marker_velocity,bad_frames_agg,params);

%% setup cluster properties
opts.whiten = 0;
opts.frameNormalize = 0;
opts.clustermethod = 'GMM';
% num = pcuse;
opts.ds = 1; % down sampling
opts.samprate = 100;
opts.params = struct;
opts.params.samplingFreq = 100;
opts.params.numPeriods=25; %distinct number of frequencies to use
opts.params.minF = 0.3; % min freq to analyze
opts.params.maxF = 20; % max freq to analyze


% for gmm model parameters
opts.pcuse = 20;
opts.numclusters = 100;
opts.lambda = 0.1; % regularization


mmm=2;

  fprintf('clipping and hipass filtering marker positions \n')
    params_clip.fps = fps;
    [markers_clipped,clipped_index_agg{mmm}] = hipass_clip(markers_preproc_aligned,badframe_inds_agg{mmm},params_clip);
    
    figure(4444)
    plot(markers_clipped.HeadF(:,3))
    
    agg_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));
    
    for ll = cluster_marker_inds      
        agg_features = cat(1,agg_features,markers_clipped.(marker_names{ll})');
    end