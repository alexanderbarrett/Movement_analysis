function descriptorstruct = get_mocap_files_table(index_in,ratname)
[num,text,mocaptable] = xlsread('Y:\Jesse\Data\Motionanalysis_captures\fileanalysis_logbook.xlsx',ratname);
mocapnames = mocaptable(:,1);
mocapnames(find(cellfun(@sum,cellfun(@isnan,mocapnames, 'UniformOutput', false)))) = [];
if index_in+1<=size(mocaptable,2)
descriptorstruct = cell2struct(mocaptable(1:numel(mocapnames),index_in+1), mocapnames, 1);
else
    fprintf('no column for that index \n')
end
end