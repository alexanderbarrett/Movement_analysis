function [agg_dprimes,agg_names] = get_feature_dprimes(ML_features,outputvector,fieldnames_beh,actualnames_beh)
%fieldnames_beh are the good fieldnames
features_to_compare = {'ja_dyadic_spectrograms',...
    'pose_score','appearance_features_agg_score_whitened','spectrogram_pcs_wl_head_angle','spectrogram_pcs_wl_trunk_angle',...
    'spectrogram_pcs_trunk_angle','spectrogram_pcs_head_angle',...
    'absolute_velocity_trunk_abs_100','absolute_std_velocity_trunk_abs_100',...
    'rel_velocity_hipR_abs_100','rel_velocity_hipL_abs_100','rel_std_velocity_hipR_abs_100','rel_std_velocity_hipL_abs_100',...
   'rel_velocity_head_abs_100','rel_std_velocity_head_abs_100','rel_velocity_trunk_abs_100','rel_std_velocity_trunk_abs_100',...
   'rel_velocity_trunk_z_100','rel_velocity_hipR_abs_33','rel_velocity_hipL_abs_33','rel_std_velocity_hipR_abs_33','rel_std_velocity_hipL_abs_33',...
   'rel_velocity_head_abs_33','rel_std_velocity_head_abs_33','rel_velocity_trunk_abs_33','rel_std_velocity_trunk_abs_33',...
   'rel_velocity_trunk_z_33','absolute_velocity_trunk_abs_33','absolute_std_velocity_trunk_abs_33',...
   'rel_velocity_hipR_abs_300','rel_velocity_hipL_abs_300','rel_std_velocity_hipR_abs_300','rel_std_velocity_hipL_abs_300',...
   'rel_velocity_head_abs_300','rel_std_velocity_head_abs_300','rel_velocity_trunk_abs_300','rel_std_velocity_trunk_abs_300',...
   'rel_velocity_trunk_z_300','absolute_velocity_trunk_abs_300','absolute_std_velocity_trunk_abs_300','rel_velocity_head_z_100','rel_std_velocity_head_z_100'};

num_unique_fieldnames = numel(fieldnames_beh);

agg_dprimes = zeros(0,num_unique_fieldnames,num_unique_fieldnames);
agg_means = zeros(0,num_unique_fieldnames);
agg_skews = zeros(0,1);

agg_wt_dprimes = zeros(0,num_unique_fieldnames,num_unique_fieldnames);

agg_names = {};

profile on
for ll = 1:numel(features_to_compare)
    feat_size = min(size(ML_features.(features_to_compare{ll})));
    for mm = 1:feat_size
        feat_diff = zeros(num_unique_fieldnames,num_unique_fieldnames);
                feat_diff_wt = zeros(num_unique_fieldnames,num_unique_fieldnames);

                feat_mean = zeros(numel(fieldnames_beh),1);

        agg_names{end+1} = strcat(features_to_compare{ll},num2str(mm));
                          overall_std =         nanstd(ML_features.(features_to_compare{ll})(:,mm),[],1);
    fprintf('comparing feature %f sub feat %f \n',ll,mm);
agg_skews = cat(1,agg_skews,skewness(ML_features.(features_to_compare{ll})(:,mm)));
        for i=reshape((fieldnames_beh),1,[])
            beh_ind_1 = find(outputvector == i);
                            mean_1 =         nanmean(ML_features.(features_to_compare{ll})(beh_ind_1,mm),1);
                            std_1 =         nanstd(ML_features.(features_to_compare{ll})(beh_ind_1,mm),[],1);

                            feat_mean(find(fieldnames_beh==i),1) = mean_1./overall_std;

            for j=reshape((fieldnames_beh),1,[])
                beh_ind_2 = find(outputvector == j);
                %% loop over features to get the differences
                
                mean_2 =         nanmean(ML_features.(features_to_compare{ll})(beh_ind_2,mm),1);
                std_2 =         nanstd(ML_features.(features_to_compare{ll})(beh_ind_2,mm),[],1);

              %  overall_std =         nanstd(ML_features.(features_to_compare{ll})(cat(2,beh_ind_1,beh_ind_2),mm),1);

                if (i ~= j)
                feat_diff(find(fieldnames_beh==i),find(fieldnames_beh==j)) = abs(mean_1./std_1-mean_2./std_2);
                                feat_diff_wt(find(fieldnames_beh==i),find(fieldnames_beh==j)) = abs(mean_1./overall_std)*abs(mean_1./std_1-mean_2./std_2);

                else
                    feat_diff(find(fieldnames_beh==i),find(fieldnames_beh==j)) = abs(mean_1)./overall_std;
                    feat_diff_wt(find(fieldnames_beh==i),find(fieldnames_beh==j)) = 0;

                end
            end
        end
        agg_dprimes = cat(1,agg_dprimes,reshape(feat_diff,1,size(feat_diff,1),size(feat_diff,2)));
        agg_wt_dprimes = cat(1,agg_wt_dprimes,reshape(feat_diff_wt,1,size(feat_diff,1),size(feat_diff,2)));
        agg_means = cat(1,agg_means,feat_mean');

    end
end
    
%% aggregate dpimes
figure(41)
imagesc(squeeze(sum(agg_wt_dprimes,1)))
    caxis([0 50])
    set(gca,'XTick',1:numel(fieldnames_beh),'XTickLabels',actualnames_beh(fieldnames_beh+1),'YTick',1:numel(fieldnames_beh),'YTickLabels',actualnames_beh(fieldnames_beh+1))
xtickangle(90)

%% indiv dprimes
plotsel = 1;
    figure(40)
    imagesc(squeeze(agg_wt_dprimes(plotsel,:,:)))
    colorbar;
    caxis([0 5])
    set(gca,'XTick',1:numel(fieldnames_beh),'XTickLabels',actualnames_beh(fieldnames_beh+1),'YTick',1:numel(fieldnames_beh),'YTickLabels',actualnames_beh(fieldnames_beh+1))
xtickangle(90)
ntitle(agg_names{plotsel})

figure(38)
auc_1 = zeros(size(agg_dprimes,1),1);
auc_2 = zeros(size(agg_dprimes,1),1);

numplot = size(agg_dprimes,1);
hold on
for mm =1:numplot%size(agg_dprimes,1)
    valshere = squeeze(agg_dprimes(mm,:,:));
    [f,x] = ecdf(valshere(find(triu(valshere,1))));
    plot(x,f);
    auc_1(mm) = sum((1-f(2:end)).*diff(x));
    restricted_ind = find(x>2);
    auc_2(mm) = sum(1-f(restricted_ind(2:end)).*diff(x(restricted_ind)));
end
legend(agg_names(1:numplot))

%% look at best performing and worst performing features
[vals,inds] = sort(auc_2,'DESCEND')
agg_names(inds(vals>0))'
%cat(2,vals(vals>0), agg_names(inds(vals>0))')

[vals,inds] = sort(auc_1,'DESCEND')
agg_names(inds)'

figure(49)
plot(vals)

%% look at best and worst discriminated pairs





