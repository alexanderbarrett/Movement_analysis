%this script takes in a list of frame intensities from an imagej z-stack
%and spits out a list of frames to remove/that are saturated. It calculates
%the frames to remove by looking at the change in intensity

%Jesse Marshall 1-16-2014


directory = '';


filename = 'm782-saturated-frames.txt';
file_save = 'on_off_behavior.txt';


delimiterIn = '\t';%'\t'
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

means = A.data(:,3);
means = means-mean(means);
differences = diff(means);

figure(1)
plot(differences)


figure(2)
plot(means)


%will need to add 1
pos_thresh = 60;
neg_thresh = 60;
on_times = find(differences > pos_thresh);
off_times = find(differences < -neg_thresh);

time_window = 100;

%off_refound = zeros(1,length(off_times));
%on_refound = zeros(1,length(off_times));

for j=1:length(off_times)
  ind_close =  find(off_times < (off_times(j) + time_window) & off_times > (off_times(j)-time_window));
    off_times(j) = max(off_times(ind_close));
end
off_times_2 = unique(off_times);


for j=1:length(on_times)
  ind_close =  find(on_times < (on_times(j) + time_window) & on_times > (on_times(j)-time_window));
    on_times(j) = min(on_times(ind_close));
end
on_times_2 = unique(on_times);

for j=1:length(on_times_2)
    if (min(abs(on_times_2(j)-off_times_2)) > 15)
   on_times_2(j) = 0;
    end
end

for j=1:length(off_times_2)
    if (min(abs(on_times_2-off_times_2(j))) > 15)
   off_times_2(j) = 0;
    end
    
    
end

on_times_2 = on_times_2(on_times_2>0);
off_times_2 = off_times_2(off_times_2>0);



remove_ind_on = [];
for j=1:length(on_times_2)
    if (min(abs(on_times_2(j)-off_times_2)) <5)
   remove_ind_on = [remove_ind_on j];
    end    
end

   remove_ind_off = [];

for j=1:length(off_times_2)
    if (min(abs(on_times_2-off_times_2(j))) <5)
   remove_ind_off = [remove_ind_off j];
    end 
end

on_times_2(remove_ind_on) = 0;
off_times_2(remove_ind_off) = 0;

on_times_2 = on_times_2(on_times_2>0);
off_times_2 = off_times_2(off_times_2>0);



times = [on_times_2 off_times_2];


A = cat(2,directory,file_save);
save(A,'times','-ascii', '-tabs')


%on_off_list



%end_time = 37259;




%csvwrite('edited_frame_list.dat',frames)
