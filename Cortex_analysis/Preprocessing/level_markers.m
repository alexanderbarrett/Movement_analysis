function markers_preproc_leveled = level_markers(markers_preproc)

fn_here = fieldnames(markers_preproc);
markers_preproc_leveled = markers_preproc;

%% visualize initial
select_mm = [5,6,7,8,9];
intercept_vals = zeros(1,numel(select_mm));
x_ind = 1;
y_ind = 2;
edge_offset = 50;
middle_offset = 20;

do_plots = 0;
default_val = 1;
for mm = select_mm
time_subset_line = setxor(1:size(markers_preproc.(fn_here{mm}),1),find(sum(markers_preproc.(fn_here{mm}),2)==0));
if (do_plots)
figure(80+mm)
subplot(1,3,1)
line(markers_preproc.(fn_here{mm})(time_subset_line,1),markers_preproc.(fn_here{mm})(time_subset_line,2))
subplot(1,3,2)
line(markers_preproc.(fn_here{mm})(time_subset_line,2),markers_preproc.(fn_here{mm})(time_subset_line,3))
subplot(1,3,3)
line(markers_preproc.(fn_here{mm})(time_subset_line,1),markers_preproc.(fn_here{mm})(time_subset_line,3))
end

x_range = intersect(find(markers_preproc.(fn_here{mm})(time_subset_line,x_ind)>min(markers_preproc.(fn_here{mm})(time_subset_line,x_ind))+edge_offset),...
    find(markers_preproc.(fn_here{mm})(time_subset_line,x_ind)<max(markers_preproc.(fn_here{mm})(time_subset_line,x_ind))-edge_offset));

y_range = intersect(find(markers_preproc.(fn_here{mm})(time_subset_line,y_ind)>min(markers_preproc.(fn_here{mm})(time_subset_line,y_ind))+edge_offset),...
    find(markers_preproc.(fn_here{mm})(time_subset_line,y_ind)<max(markers_preproc.(fn_here{mm})(time_subset_line,y_ind))-edge_offset));
if (numel(x_range) && numel(y_range))

x_model = fitlm(markers_preproc.(fn_here{mm})(time_subset_line(x_range),x_ind),markers_preproc.(fn_here{mm})(time_subset_line(x_range),3));
y_model = fitlm(markers_preproc.(fn_here{mm})(time_subset_line(y_range),y_ind),markers_preproc.(fn_here{mm})(time_subset_line(y_range),3));

intercept_vals(1,find(select_mm == mm)) = x_model.Coefficients.Estimate(2);
intercept_vals(2,find(select_mm == mm)) = y_model.Coefficients.Estimate(2);
else
  
intercept_vals(1,find(select_mm == mm)) = nan;  
intercept_vals(2,find(select_mm == mm)) = nan;  

end
end

fprintf('first round of intercepts');
intercept_vals

%get the rotation about the y_axis
y_axis_ang = atan(mean(intercept_vals(1,:),2));
if ~isnan(y_axis_ang)
rotation_y = zeros(3,3);
rotation_y(1,:) = [cos(y_axis_ang) 0 -sin(y_axis_ang)];
rotation_y(2,:) = [0 1 0];
rotation_y(3,:) = [sin(y_axis_ang) 0 cos(y_axis_ang)];
%rotated_temp = mtimesx(markers_preproc.(fn_here{mm}),rotation_y);

for mm = 1:numel(fn_here)
markers_preproc_leveled.(fn_here{mm}) = mtimesx(markers_preproc_leveled.(fn_here{mm}),rotation_y);
end
end
% get the rotation about the x-axis
for mm = select_mm
time_subset_line = setxor(1:size(markers_preproc.(fn_here{mm}),1),find(sum(markers_preproc.(fn_here{mm}),2)==0));

x_range = intersect(find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind)>min(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind))+edge_offset),...
    find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind)<max(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind))-edge_offset));

y_range = intersect(find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind)>min(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind))+edge_offset),...
    find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind)<max(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind))-edge_offset));
if (numel(x_range) && numel(y_range))

x_model = fitlm(markers_preproc_leveled.(fn_here{mm})(time_subset_line(x_range),x_ind),markers_preproc_leveled.(fn_here{mm})(time_subset_line(x_range),3));
y_model = fitlm(markers_preproc_leveled.(fn_here{mm})(time_subset_line(y_range),y_ind),markers_preproc_leveled.(fn_here{mm})(time_subset_line(y_range),3));

intercept_vals(1,find(select_mm == mm)) = x_model.Coefficients.Estimate(2);
intercept_vals(2,find(select_mm == mm)) = y_model.Coefficients.Estimate(2);
else
    intercept_vals(1,find(select_mm == mm)) = nan;  
intercept_vals(2,find(select_mm == mm)) = nan;  
end
end

fprintf('secound round of intercepts');

intercept_vals

x_axis_ang = -atan(mean(intercept_vals(2,:),2));
x_axis_ang
if ~isnan(x_axis_ang)

rotation_x = zeros(3,3);
rotation_x(1,:) = [1 0 0];
rotation_x(2,:) =  [0 cos(x_axis_ang) sin(x_axis_ang)];
rotation_x(3,:) = [0 -sin(x_axis_ang) cos(x_axis_ang)];

for mm = 1:numel(fn_here)
markers_preproc_leveled.(fn_here{mm}) = mtimesx(markers_preproc_leveled.(fn_here{mm}),rotation_x);
end
end

if (do_plots)
figure(90)
subplot(1,3,1)
line(markers_preproc_leveled.(fn_here{mm})(time_subset_line,1),markers_preproc_leveled.(fn_here{mm})(time_subset_line,2))
subplot(1,3,2)
line(markers_preproc_leveled.(fn_here{mm})(time_subset_line,2),markers_preproc_leveled.(fn_here{mm})(time_subset_line,3))
subplot(1,3,3)
line(markers_preproc_leveled.(fn_here{mm})(time_subset_line,1),markers_preproc_leveled.(fn_here{mm})(time_subset_line,3))
end



% get the rotation about the x-axis
for mm = select_mm
   
time_subset_line = setxor(1:size(markers_preproc.(fn_here{mm}),1),find(sum(markers_preproc.(fn_here{mm}),2)==0));

x_range = intersect(find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind)>min(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind))+edge_offset),...
    find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind)<max(markers_preproc_leveled.(fn_here{mm})(time_subset_line,x_ind))-edge_offset));

y_range = intersect(find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind)>min(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind))+edge_offset),...
    find(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind)<max(markers_preproc_leveled.(fn_here{mm})(time_subset_line,y_ind))-edge_offset));
if (numel(x_range) && numel(y_range))

x_model = fitlm(markers_preproc_leveled.(fn_here{mm})(time_subset_line(x_range),x_ind),markers_preproc_leveled.(fn_here{mm})(time_subset_line(x_range),3));
y_model = fitlm(markers_preproc_leveled.(fn_here{mm})(time_subset_line(y_range),y_ind),markers_preproc_leveled.(fn_here{mm})(time_subset_line(y_range),3));

intercept_vals(1,find(select_mm == mm)) = x_model.Coefficients.Estimate(2);
intercept_vals(2,find(select_mm == mm)) = y_model.Coefficients.Estimate(2);
else
    intercept_vals(1,find(select_mm == mm)) = nan;  
intercept_vals(2,find(select_mm == mm)) = nan;  
end
end
fprintf('final intercepts \n')
intercept_vals



%% get initial slopes -- exclude ends of the arena where animal rears more
