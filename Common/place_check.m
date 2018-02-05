function vargout = place_check(picture,centroids,traces)


dayday =3;
centroids = experiment_data{1}.centroid_good{dayday};
portdist = experiment_data{1}.port_dist{dayday};
trace_length = experiment_data{1}.length_trace{dayday};
%traces = experiment_data{1}.length_trace{dayday}

centroid_x = conv(centroids(1,:),ones(1,4)./4,'same');
centroid_y = conv(centroids(2,:),ones(1,4)./4,'same');

tracking_length = min(4*trace_length,numel(centroid_x));
tracking_length_5hz = tracking_length./4;

centroid_x = centroid_x(1:4:tracking_length);
centroid_y = centroid_y(1:4:tracking_length);
portdist = portdist(1:4:tracking_length);

fprintf(experiment_data{1}.name_tag{dayday})


figure(1)
hold on
line(centroid_y,centroid_x)
hold off

max_y = max(centroid_y);
max_x = max(centroid_x);
num_bins = 25;
x_edges = 0:max_x/num_bins:max_x;
y_edges = 0:max_y/num_bins:max_y;
figure(999)
histmat = hist2(centroid_y,centroid_x,y_edges,x_edges);
imagesc(y_edges,x_edges,histmat'); colorbar ; axis square tight ;

sum_trace = zeros(1,trace_length);


for cellid = experiment_data{1}.day_id{dayday}
    sum_trace = sum_trace + experiment_data{cellid}.events{dayday}';
end
sum_trace = sum_trace(1:trace_length);

value_hist = zeros(num_bins,num_bins);
%there has to be a meshgrid way
for xk = 1:num_bins
    for yk=1:num_bins
   place_times = intersect(find((y_edges(yk) < centroid_y )& (centroid_y <  y_edges(yk+1))),  find((x_edges(xk) < centroid_x) & (centroid_x <  x_edges(xk+1) )))   ;
   if (place_times)
      value_hist(yk,xk) = mean(sum_trace(place_times)); 
   end
   
    end
end

%look at center of mass weighted 
cell_pds = zeros(1,numel(experiment_data{1}.day_id{dayday}));
cell_stds = zeros(1,numel(experiment_data{1}.day_id{dayday}));

for cellid = experiment_data{1}.day_id{dayday}
cell_times = find(experiment_data{cellid}.events{dayday} == 1);
cell_times(cell_times>tracking_length_5hz) = tracking_length_5hz;
std_y = std(centroid_y(cell_times));
std_x = std(centroid_x(cell_times));
cellid
if (sum(experiment_data{cellid}.events{dayday}) > 10)
cell_pds(find(experiment_data{1}.day_id{dayday}) == cellid) = mean(portdist(cell_times));
cell_stds(find(experiment_data{1}.day_id{dayday}) == cellid) = std_y+std_x;
end
end

cell_stds(cell_stds == 0) = 10000000;
cell_pds(cell_pds == 0) = 10000000;
[val,ind] = sort(cell_stds,'Ascend');


for cellid = ind(1:50)
    number = find(ind(1:50) == cellid);
    cell_times = find(experiment_data{cellid}.events{dayday} == 1);
cell_times(cell_times>tracking_length_5hz) = tracking_length_5hz;
%plot top
figure(1000+number)
subplot(2,1,1)
line(centroid_y,centroid_x)
hold on
plot(centroid_y(cell_times),centroid_x(cell_times),'r+')
hold off

subplot(2,1,2)
plot(2.*experiment_data{ind(cellid)}.smoothed_events{dayday});
hold on
plot(experiment_data{1}.lick_trace{dayday},'r');
plot(2.*conv(experiment_data{1}.CS_trace{dayday},ones(1,10),'same'),'k');
hold off
end


cellinspect = 19;
figure(2000)



%plot an individual IC
figure(800)
imagesc(y_edges,x_edges,value_hist')

 result = 100;
    while (result > 0)
        
prompt = 'Enter IC number to Examine [0 for exit] ';
result = input(prompt);
        if (result >0)

cell_times = find(experiment_data{result}.events{dayday} == 1);
cell_times(cell_times>tracking_length) = 0;
figure(1000)
line(centroid_y,centroid_x)
hold on
plot(centroid_y(cell_times),centroid_x(cell_times),'r+')
hold off


        else 
            return
        end
    end
end