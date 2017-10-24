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
title('Pre')



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
title('Post')


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



end