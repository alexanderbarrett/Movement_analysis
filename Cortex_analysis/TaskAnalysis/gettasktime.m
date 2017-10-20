function [licktime,levertime] = gettasktime(mocapstruct)

levertime = zeros(1,numel(mocapstruct.lever_thresholded));
temp = conv(mocapstruct.lever_thresholded,ones(1,2*mocapstruct.fps),'same');
levertime(temp>0) = 1;


licktime = zeros(1,numel(mocapstruct.resample_analog(:,5)));
temp = conv(mocapstruct.resample_analog(:,5),ones(1,2*mocapstruct.fps),'same');
licktime(temp>0) = 1;


end