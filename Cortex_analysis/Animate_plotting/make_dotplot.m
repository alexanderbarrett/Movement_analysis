function make_dotplot(annotation_vec,fieldnames_beh,ML_features,beh_subset)
feature_list = {'LGroom','high_rear','rel_velocity_head_abs_100','rel_velocity_head_abs_33',...
    'absolute_velocity_trunk_abs_100','rel_velocity_hipL_abs_100','rel_velocity_hipR_abs_100','rel_velocity_trunk_abs_100'};
feature_names = {'LGroom','high rear','Head Vel','Head Velocity (faster)',...
    'Trunk Velocity','L Hip Velocity','R Hip Velocity','Rel. Trunk Velocity'};

% get the max val and std of the chosen features for z scoring
feature_means = zeros(1,numel(feature_list)); 
feature_stds = zeros(1,numel(feature_list)); 

for ll = 1:numel(feature_list)
    if numel(size(ML_features.(feature_list{ll}))==2) && ...
            size(ML_features.(feature_list{ll}),1)>size(ML_features.(feature_list{ll}),2)
    ML_features.(feature_list{ll}) = ML_features.(feature_list{ll})';
    end
      feature_means(ll) = nanmean(ML_features.(feature_list{ll})(1,:));
            feature_stds(ll) = nanstd(ML_features.(feature_list{ll})(1,:));

  end
% loop over all behaviors and maek the plot
behavior_ids = 0:max(annotation_vec);%unique(annotation_vec);

behavior_ids(isnan(behavior_ids)) = [];
num_beh = numel(behavior_ids);
cmap = jet(100);

beh_val = zeros(num_beh,numel(feature_list));
beh_std = zeros(num_beh,numel(feature_list));
beh_color = zeros(num_beh,numel(feature_list));


if nargin == 4
    num_beh = numel(beh_subset);
else
    beh_subset = 1:num_beh;
end
beh_subset = reshape(beh_subset,1,[]);

axisnames = cell(1,0);

for mm = beh_subset
    size(axisnames,2);
      axisnames{1,size(axisnames,2)+1} = fieldnames_beh{(mm)+1};
    behaviors_inds = find(annotation_vec == (mm));
    for ll = 1:numel(feature_list)
        
beh_val(mm,ll) = (nanmean(ML_features.(feature_list{ll})(1,behaviors_inds))-feature_means(ll))./feature_stds(ll);
beh_std(mm,ll) = (nanstd(ML_features.(feature_list{ll})(1,behaviors_inds))-feature_stds(ll))./feature_stds(ll);
    end
end

%% normalize the values
beh_val(beh_val>2) = 2;
beh_val(beh_val<-2) = -2;
beh_val = beh_val+2;
scaled_beh_val = beh_val;
scaled_beh_std = beh_std;

for ll = 1:numel(feature_list)
scaled_beh_val(:,ll) = 1+(beh_val(:,ll)-min(squeeze(beh_val(:,ll))))*9./...
    (max(squeeze(beh_val(:,ll)))-min((squeeze(beh_val(:,ll)))));

scaled_beh_std(:,ll) = 1+(beh_std(:,ll)-min(squeeze(beh_std(:,ll))))*9./...
    (max(squeeze(beh_std(:,ll)))-min((squeeze(beh_std(:,ll)))));
end

% beh_color = beh_std+1;
% beh_color(beh_color<0) = 0;
% beh_color(beh_color>2) = 2;
% beh_color = round(beh_color*50);
beh_color = round(scaled_beh_std*10);



    figure (111)

for mm = beh_subset
    for ll = 1:numel(feature_list)
        colorhere = 'k';%cmap(beh_color(mm,ll),:);
        if ~isnan(scaled_beh_val(mm,ll))
        plot(ll,find(beh_subset==mm),'o','MarkerSize',(scaled_beh_val(mm,ll)),'MarkerFaceColor',colorhere,'MarkerEdgeColor',colorhere)
        end
%plot(ll,mm,'o','MarkerSize',beh_mult(ll)*(beh_val(mm,ll))./2+1,'MarkerFaceColor','k','MarkerEdgeColor','k')
%,cmap(beh_color+1,:)
%cmap(beh_color+1,:)
hold on
end
    end

xlim([0 9])
set(gca,'XTick',[1:numel(feature_list)],'XTickLabels',feature_names,'YTick',[1:num_beh],'YTickLabels',axisnames)
colormap(jet)
colorbar;

% %declare colormaps
% num_beh = numel(fieldnames_beh);
% ethocolors = hsv(num_beh);
% %get behavior to number mapping
%    frames_beh = find(annotation_vec == mm); 
%       if numel(frames_beh)
%    figure(67)
%   h= bar( frames_beh,ones(1,numel( frames_beh)),1)
%   set(h,'FaceColor',ethocolors(mm+1,:),'EdgeColor','none')
%    hold on
%    legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};
%    end  
% 
% legend(legendnames)
% 

%pose ethogram


% head ethogram

% 