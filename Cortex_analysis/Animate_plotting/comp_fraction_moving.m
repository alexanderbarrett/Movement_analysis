function comp_fraction_moving(mocapstruct)
frac_move = numel(mocapstruct.move_frames)./(numel(mocapstruct.move_frames)+numel(mocapstruct.rest_frames));
fprintf('fraction of time moving: %f \n',frac_move);
end