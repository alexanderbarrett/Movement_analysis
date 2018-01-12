function [descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory)

mocapfilestruct = loadmocapfilestruct(descriptor_struct.Rat,mocapmasterdirectory);
% 
% switch descriptor_struct.Rat
%     case 'Vicon3'
% mocapfilestruct = loadmocapfilestruct('Vicon3',mocapmasterdirectory);
% %mocapfilestruct.mocapdir = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\';
% 
%     case 'Vicon8'
% mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);
%     case 'JDM25'
% mocapfilestruct = loadmocapfilestruct('JDM25',mocapmasterdirectory);
%    case 'JDM32'
% mocapfilestruct = loadmocapfilestruct('JDM32',mocapmasterdirectory);
% case 'JDM33'
% mocapfilestruct = loadmocapfilestruct('JDM33',mocapmasterdirectory);
% case 'JDM27'
% mocapfilestruct = loadmocapfilestruct('JDM33',mocapmasterdirectory);
%   case 'JDM21'
% mocapfilestruct = loadmocapfilestruct('JDM21',mocapmasterdirectory);
% end
mocapvideodirectory = [];
%% get the desired files


mocapfilearray=[];
 mocapfiletimes = [];
 if (ischar(descriptor_struct.Days))
 days_search = str2num(descriptor_struct.Days);
 else 
      days_search = (descriptor_struct.Days);
 end
 for day_find =  days_search
    % day_find = str2num(day_find);
descriptor_struct_1.day = (day_find);
if numel(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_find})>1
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_find},descriptor_struct.Searchtag)));
else
 good_inds = find(numel(strfind(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_find},descriptor_struct.Searchtag)));
end
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct.Condition).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct.Condition).threshcrossings{day_find}(good_inds);
switch descriptor_struct.Task
    case 0
nontaskfiles = find(mocapthreshcrossings<descriptor_struct.TaskThreshold);
case 1
nontaskfiles = find(mocapthreshcrossings>descriptor_struct.TaskThreshold);
case 2
nontaskfiles = good_inds;
end
good_inds = intersect(good_inds,nontaskfiles);

if numel(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_find}(good_inds)')
mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct.Condition).mocapdatecreate{day_find}(good_inds)');
else 
    fprintf('no files found \n')
end



 end
 
 if (descriptor_struct.HourRange_early~=0)
     if (descriptor_struct.HourRange_early > descriptor_struct.HourRange_late)
  goodinds = union(find(datetime(mocapfiletimes).Hour>descriptor_struct.HourRange_early) , find(datetime(mocapfiletimes).Hour<descriptor_struct.HourRange_late));
     else
         % during the day
    goodinds = intersect(find(datetime(mocapfiletimes).Hour>descriptor_struct.HourRange_early) , find(datetime(mocapfiletimes).Hour<descriptor_struct.HourRange_late));      
     end
 else         
     goodinds = 1:numel(mocapfiletimes);
 end
  numsample = min(numel(goodinds),descriptor_struct.NumFiles);
 goodinds_sample = goodinds(1:numsample);%randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
 




end