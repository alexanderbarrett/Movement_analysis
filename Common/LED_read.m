% filename = 'on_off_scope.txt';
% filename_2 = '130209-F506-PAV2-TTL.txt';
% file_out = 'edited_frame_list_scope.dat';

filename = 'on_off_behavior.txt';
%filename_2 = '130210-F525-PAV1-TTL-TIMES.txt';
file_out = 'edited_frame_list_behavior.dat';
filename_2 = 'm782-saturated-frames.txt';


delimiterIn = '\t';%'\t'
headerlinesIn = 1;


directory = '';

delimiterIn_2 = '\t';%'\t'
headerlinesIn_2 = 1;
A_2 = importdata(filename_2,delimiterIn_2,headerlinesIn_2);
means = A_2.data(:,1);

A = importdata(filename,delimiterIn,headerlinesIn);

trial_index = 1:25;
on_times = A.data(:,1);
off_times = A.data(:,2);

ind = [];
for i = 1:length(on_times)
   ind = cat(2,ind,(on_times(i):off_times(i)));
end

end_time = max(means);

frames = 1:end_time;
frames = setxor(frames,ind);


csvwrite(cat(2,directory,file_out),frames)
