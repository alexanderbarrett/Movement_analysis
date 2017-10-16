% function to copy marker files

marker_file_1 = 'E:\Bence\Data\Motionanalysis_captures\JDM25\20171010\marker_old.csv';
marker_file_2 = 'E:\Bence\Data\Motionanalysis_captures\JDM25\20171010\marker_new.csv';


links_to_reset = [6,14,17,20];
T = readtable(marker_file_1);
TT = table2array(T);
TT_new = TT;
rows_to_mod = [];
for mm = 1:size(TT,1)
    if ismember(TT(mm,3),links_to_reset) || ismember(TT(mm,4),links_to_reset)
rows_to_mod = cat(1,rows_to_mod,mm);
mean_val = mean([TT(mm,5) TT(mm,6)]);
TT_new(mm,5:6) = mean_val;
    end
end

csvwrite(marker_file_2,TT_new);