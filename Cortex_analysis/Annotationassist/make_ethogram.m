function make_ethogram(annotation_vec,fieldnames_beh)
%declare colormaps
num_beh = numel(fieldnames_beh);
ethocolors = hsv(num_beh);
%get behavior to number mapping
legendnames = cell(1,0);
for mm = 0:num_beh-1
   frames_beh = find(annotation_vec == mm); 
      if numel(frames_beh)
   figure(67)
  h= bar( frames_beh,ones(1,numel( frames_beh)),1)
  set(h,'FaceColor',ethocolors(mm+1,:),'EdgeColor','none')
   hold on
   legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};
   end  
end
legend(legendnames)


%pose ethogram


% head ethogram

% 