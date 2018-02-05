%convert from spikeE
trace_ind = 3;
image_ind = 3;
temp1 = struct2cell(SpikeTraceData);
temp2 = temp1(trace_ind,1,:);
temp3 = squeeze(cell2mat(temp2));
IcaTraces = permute(temp3,[2 1]);
save('ICtraces_sorted.mat','IcaTraces')


temp1 = struct2cell(SpikeImageData);
temp2 = temp1(image_ind,1,:);
temp3 = squeeze(cell2mat(temp2));
IcaFilters = permute(temp3,[3 1 2]);
save('ICfilters_sorted','IcaFilters')