function make_dotplot(annotation_vec,fieldnames_beh,ML_features)
feature_list = {'LGroom','RGroom','high_rear','rel_velocity_head_abs_100',...
    'absolute_velocity_trunk_abs_100','rel_velocity_hipL_abs_100','rel_velocity_hipR_abs_100','rel_velocity_trunk_abs_100'};


% get the max val and std of the chosen features for z scoring
feature_means = zeros(1,numel(feature_list)); 
feature_stds = zeros(1,numel(feature_list)); 
axisnames = cell(1,0);


beh_mult = 4*ones(1,numel(feature_list));
beh_mult(4:end) = 8;

for ll = 1:numel(feature_list)
    if numel(size(ML_features.(feature_list{ll}))==2) && ...
            size(ML_features.(feature_list{ll}),1)>size(ML_features.(feature_list{ll}),2)
    ML_features.(feature_list{ll}) = ML_features.(feature_list{ll})';
    end
      feature_means(ll) = nanmean(ML_features.(feature_list{ll})(1,:));
            feature_stds(ll) = nanstd(ML_features.(feature_list{ll})(1,:));

  end
% loop over all behaviors and maek the plot
behavior_ids = unique(annotation_vec);
behavior_ids(isnan(behavior_ids)) = [];
num_beh = numel(behavior_ids);
cmap = jet(100);


for mm = 1:num_beh
      axisnames{1,size(axisnames,2)+1} = fieldnames_beh{behavior_ids(mm)+1};
    behaviors_inds = find(annotation_vec == behavior_ids(mm));
    for ll = 1:numel(feature_list)
        
beh_val = (nanmean(ML_features.(feature_list{ll})(1,behaviors_inds))-feature_means(ll))./feature_stds(ll);
beh_std = (nanstd(ML_features.(feature_list{ll})(1,behaviors_inds))-feature_stds(ll))./feature_stds(ll);

beh_val(beh_val>2) = 2;
beh_val(beh_val<-2) = -2;
beh_val = beh_val+2;
%normalize the color of the marker
beh_color = beh_std+1;
beh_color(beh_color<0) = 0;
beh_color(beh_color>2) = 2;

beh_color = round(beh_color*50);

    figure (111)
plot(ll,mm,'o','MarkerSize',beh_mult(ll)*(beh_val)./2+1,'MarkerFaceColor','k','MarkerEdgeColor','k')
%,cmap(beh_color+1,:)
%cmap(beh_color+1,:)
hold on
    end


end
xlim([0 9])
set(gca,'XTick',[1:numel(feature_list)],'XTickLabels',feature_list,'YTick',[1:num_beh],'YTickLabels',axisnames)
colormap(jet)
colorbar;

%declare colormaps
num_beh = numel(fieldnames_beh);
ethocolors = hsv(num_beh);
%get behavior to number mapping
   frames_beh = find(annotation_vec == mm); 
      if numel(frames_beh)
   figure(67)
  h= bar( frames_beh,ones(1,numel( frames_beh)),1)
  set(h,'FaceColor',ethocolors(mm+1,:),'EdgeColor','none')
   hold on
   legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};
   end  

legend(legendnames)


%pose ethogram


% head ethogram

% 