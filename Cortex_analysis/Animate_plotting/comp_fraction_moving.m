function comp_fraction_moving(mocapstruct)
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

xval = 100:100:10000;

figure(33)
[n,x] = hist(move_lengths,xval);
loglog(xval,n./sum(n),'r','linewidth',2)
hold on
[n2,x2] = hist(rest_lengths,xval);
loglog(xval,n2./sum(n2),'b','linewidth',2)
hold off
legend('move','rest')
box off
xlabel('time(frames)')
ylabel('fraction of bouts')


end