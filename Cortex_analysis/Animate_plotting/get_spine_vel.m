function velsum_pre_move= get_spine_vel(mocapstruct)

velsum_post = [];
difforder = 10;

for mm =1:3
    smooth_vel = medfilt1(conv(squeeze(mocapstruct.aligned_mean_position(:,mm)),ones(1,10)./10,'same'),1);
    if mm == 1
    good_ind = find(smooth_vel ~=0);
    velsum_pre =diff(smooth_vel(good_ind(1:difforder:end))).^2;
    
    good_ind_move = intersect(good_ind,mocapstruct.move_frames);
     good_ind_rest = intersect(good_ind,mocapstruct.rest_frames);
       velsum_pre_rest =diff(smooth_vel(good_ind_rest(1:difforder:end))).^2;
    velsum_pre_move =diff(smooth_vel(good_ind_move(1:difforder:end))).^2;

    else
velsum_pre = velsum_pre+diff(smooth_vel(good_ind(1:difforder:end))).^2;
velsum_pre_rest = velsum_pre_rest+diff(smooth_vel(good_ind_rest(1:difforder:end))).^2;
velsum_pre_move = velsum_pre_move+diff(smooth_vel(good_ind_move(1:difforder:end))).^2;

    end
end
velsum_pre = (mocapstruct.fps./difforder)*sqrt(velsum_pre);
velsum_pre_rest = (mocapstruct.fps./difforder)*sqrt(velsum_pre_rest);
velsum_pre_move = (mocapstruct.fps./difforder)*sqrt(velsum_pre_move);
end