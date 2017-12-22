function descriptorstruct = get_mocap_files_table(index_in)
[num,text,mocaptable] = xlsread('Y:\Jesse\Data\Motionanalysis_captures\fileanalysis_logbook.xlsx');
mocapnames = mocaptable(:,1);
descriptorstruct = cell2struct(mocaptable(:,index_in+1), mocapnames, 1);
end