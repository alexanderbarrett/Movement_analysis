function compare_comp_fraction_moving(mocapstruct, mocapstruct_post)
frac_move = numel(mocapstruct.move_frames)./(numel(mocapstruct.move_frames)+numel(mocapstruct.rest_frames));
fprintf('fraction of time moving: %f \n',frac_move);

frac_move_post = numel(mocapstruct_post.move_frames)./(numel(mocapstruct_post.move_frames)+numel(mocapstruct_post.rest_frames));
fprintf('fraction of time moving: %f \n',frac_move_post);

%% look at move
inst_label = ones(1,numel(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(:,3)));
inst_label(mocapstruct.rest_frames) = 0;
pixellist = bwconncomp(inst_label);
move_lengths = cellfun(@numel,pixellist.PixelIdxList);

%% do for rest
inst_label = ones(1,numel(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(:,3)));
inst_label(mocapstruct.move_frames) = 0;
pixellist = bwconncomp(inst_label);
rest_lengths = cellfun(@numel,pixellist.PixelIdxList);

%xval = 100:100:10000;
xval =logspace(1,5,50)
figure(33)
[n,x] = hist(move_lengths,xval);
loglog(xval./mocapstruct.fps,n./sum(n),'r','linewidth',2)
hold on
[n2,x2] = hist(rest_lengths,xval);
loglog(xval./mocapstruct.fps,n2./sum(n2),'b','linewidth',2)
hold off
legend('move','rest')
box off
xlabel('time(s)')
ylabel('fraction of bouts')
title('Pre move vs rest')



%% look at move
inst_label = ones(1,numel(mocapstruct_post.markers_preproc.(mocapstruct_post.markernames{1})(:,3)));
inst_label(mocapstruct_post.rest_frames) = 0;
pixellist = bwconncomp(inst_label);
move_lengths_post = cellfun(@numel,pixellist.PixelIdxList);

%% do for rest
inst_label = ones(1,numel(mocapstruct_post.markers_preproc.(mocapstruct_post.markernames{1})(:,3)));
inst_label(mocapstruct_post.move_frames) = 0;
pixellist = bwconncomp(inst_label);
rest_lengths_post = cellfun(@numel,pixellist.PixelIdxList);

%xval = 100:100:10000;
xval =logspace(1,5,50)
figure(34)
[n_post,x_post] = hist(move_lengths_post,xval);
loglog(xval./mocapstruct_post.fps,n_post./sum(n_post),'r','linewidth',2)
hold on
[n2_post,x2_post] = hist(rest_lengths_post,xval);
loglog(xval./mocapstruct_post.fps,n2_post./sum(n2_post),'b','linewidth',2)
hold off
legend('move','rest')
box off
xlabel('time(s)')
ylabel('fraction of bouts')
title('Post move vs rest')


%% compare the two
figure(35)
subplot(2,1,1)
loglog(xval./mocapstruct_post.fps,n./sum(n),'r','linewidth',2)
hold on
loglog(xval./mocapstruct_post.fps,n_post./sum(n_post),'b','linewidth',2)
hold off
legend('Pre','Post')
box off
xlabel('time(s)')
ylabel('fraction of bouts')
title('Move')


subplot(2,1,2)
loglog(xval./mocapstruct_post.fps,n2./sum(n2),'r','linewidth',2)
hold on
loglog(xval./mocapstruct_post.fps,n2_post./sum(n2_post),'b','linewidth',2)
hold off
legend('Pre','Post')
box off
xlabel('time(s)')
ylabel('fraction of bouts')
title('Rest')


%% look at the spine center velocity and compare across cond
velsum_post = [];
difforder = 105;

for mm =1:3
    smooth_vel = medfilt1(conv(squeeze(mocapstruct.aligned_mean_position(:,mm)),ones(1,10)./10,'same'),1);
     smooth_vel_post = medfilt1(conv(squeeze(mocapstruct_post.aligned_mean_position(:,mm)),ones(1,10)./10,'same'),1);
    if mm == 1
    good_ind = find(smooth_vel ~=0);
    velsum_pre =diff(smooth_vel(good_ind(1:difforder:end))).^2;
    
       good_ind_post = find(smooth_vel_post ~=0);
    velsum_post =diff(smooth_vel_post(good_ind_post(1:difforder:end))).^2;
    else
velsum_pre = velsum_pre+diff(smooth_vel(good_ind(1:difforder:end))).^2;
velsum_post = velsum_post+diff(smooth_vel_post(good_ind_post(1:difforder:end))).^2;
    end
end
velsum_pre = (mocapstruct_post.fps./difforder)*sqrt(velsum_pre);
velsum_post = (mocapstruct_post.fps./difforder)*sqrt(velsum_post);
xval =logspace(-1,3,100)
figure(36)
%[n_post,x_post] = hist(move_lengths_post,xval);
%loglog(xval./mocapstruct_post.fps,n_post./sum(n_post),'r','linewidth',2)
[n2_pre,x2_pre] = hist(velsum_pre,xval);
ind_use =1:numel(x2_pre); %cancel out the zero
semilogx(xval(ind_use),n2_pre(ind_use)./sum(n2_pre(ind_use)),'b','linewidth',2)
hold on
[n2_post,x2_post] = hist(velsum_post,xval);
semilogx(xval(ind_use),n2_post(ind_use)./sum(n2_post(ind_use)),'r','linewidth',2)
hold off
legend('pre','post')
box off
xlabel('velocity (mm/s)')
ylabel('fraction of time')
title('Post move vs rest')


end