function compare_comp_fraction_moving(mocapstruct, colorin,iter)
frac_move = numel(mocapstruct.move_frames)./(numel(mocapstruct.move_frames)+numel(mocapstruct.rest_frames));
fprintf('fraction of time moving: %f \n',frac_move);

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
xval =logspace(2,5,20);
figure(33+iter)
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



%% compare the two
figure(105)
subplot(2,1,1)
hold on
loglog(xval./mocapstruct.fps,n./sum(n),'Color',colorin,'linewidth',2)
hold off
legend('Pre','Post')
box off
xlabel('Bout Length(s)')
ylabel('fraction of bouts')
set(gca,'xscale','log','yscale','log')
title('Move')


subplot(2,1,2)
hold on
loglog(xval./mocapstruct.fps,n2./sum(n2),'Color',colorin,'linewidth',2)
set(gca,'xscale','log','yscale','log')
legend('Pre','Post')
box off
xlabel('Bout Length(s)')
ylabel('fraction of bouts')
title('Rest')


%% look at the spine center velocity and compare across cond
velsum_post = [];
difforder = 20;
 badframes = [];
for mm =1:3
    badframes =  (cat(1,badframes,find(mocapstruct.aligned_mean_position(:,mm) == 0)));
    badframes = (cat(1,badframes,find(diff(mocapstruct.aligned_mean_position(:,mm),1)>10)));
end
badframes=unique(bsxfun(@plus,badframes,-2*difforder:2*difforder));

induse = [1:difforder:size(mocapstruct.aligned_mean_position,1)];
[~,badind,~] = intersect(induse,badframes);
[~,subindmove,~]= intersect(induse,mocapstruct.move_frames);
[~,subindrest,~]= intersect(induse,mocapstruct.rest_frames);


for mm =1:3
    smooth_vel =squeeze(mocapstruct.aligned_mean_position(:,mm));%(conv(squeeze(mocapstruct.aligned_mean_position(:,mm)),ones(1,10)./10,'same');
  
    sqvel = cat(1,diff(smooth_vel((induse))).^2,0);
        sqvel(badind) = nan;

    if mm == 1
    %good_ind = find(smooth_vel ~=0);
    velsum_pre =sqvel;
    velsum_pre_rest = sqvel(subindrest);
        velsum_pre_move = sqvel(subindmove);

  %  good_ind_move = intersect(good_ind,mocapstruct.move_frames);
   %  good_ind_rest = intersect(good_ind,mocapstruct.rest_frames);
      % velsum_pre_rest =diff(smooth_vel(good_ind_rest(1:difforder:end))).^2;
    %velsum_pre_move =diff(smooth_vel(good_ind_move(1:difforder:end))).^2;

    else
velsum_pre = velsum_pre+sqvel;
velsum_pre_rest = velsum_pre_rest+sqvel(subindrest);
velsum_pre_move = velsum_pre_move+sqvel(subindmove);

    end
end
velsum_pre = (mocapstruct.fps./difforder)*sqrt(velsum_pre);
velsum_pre_rest = (mocapstruct.fps./difforder)*sqrt(velsum_pre_rest);
velsum_pre_move = (mocapstruct.fps./difforder)*sqrt(velsum_pre_move);

% figure(66)
% plot(velsum_pre,'b')
% hold on
% plot(0:1./difforder:(numel(smooth_vel)-1)./difforder,smooth_vel,'r')
% hold off
num_bad = numel(badind);
xval =logspace(-1,3,100);
figure(306)
hold on
[n2_pre,x2_pre] = hist(velsum_pre,xval);
ind_use =1:numel(x2_pre); %cancel out the zero
n2_pre(2) = n2_pre(2)+num_bad;
total_frames = sum(n2_pre(ind_use))+num_bad;
semilogx(xval(ind_use),n2_pre(ind_use)./sum(n2_pre(ind_use)),'Color',colorin,'linewidth',2)
box off
set(gca,'xscale','log')
xlabel('Spine Velocity (mm/s)')
ylabel('fraction of time')
title('Spine Velocity')


figure(307)
subplot(1,2,1)
hold on
[n2_pre,x2_pre] = hist(velsum_pre_rest,xval);
ind_use =1:numel(x2_pre); %cancel out the zero
semilogx(xval(ind_use),n2_pre(ind_use)./sum(n2_pre(ind_use)),'Color',colorin,'linewidth',2)
box off
set(gca,'xscale','log')
xlabel('Spine Velocity (mm/s)')
ylabel('fraction of time')
title('Rest')

subplot(1,2,2)
hold on
[n2_pre,x2_pre] = hist(velsum_pre_move,xval);
ind_use =1:numel(x2_pre); %cancel out the zero
semilogx(xval(ind_use),n2_pre(ind_use)./sum(n2_pre(ind_use)),'Color',colorin,'linewidth',2)
box off
set(gca,'xscale','log')
xlabel('Spine Velocity (mm/s)')
ylabel('fraction of time')
title('Move')



end