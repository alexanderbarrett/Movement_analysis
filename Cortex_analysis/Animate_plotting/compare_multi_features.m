function compare_multi_features(mocapstructcellarray,markersubsets,names,comparison_type)

annotation_fields = fieldnames(mocapstructcellarray{1}.annotated_pose_struct);
comparefields = fieldnames(mocapstructcellarray{1}.annotated_pose_struct);
numbehaviors = numel(comparefields);

output_struct = cell(numel(mocapstructcellarray),numel(mocapstructcellarray),numbehaviors);
kfoldlosses = zeros(numel(mocapstructcellarray),numel(mocapstructcellarray),numbehaviors);

subcluster_names = mocapstructcellarray{1}.modular_cluster_properties.subcluster_names;

fieldsubset = [1,3:6];
for mm =1:numel(mocapstructcellarray)
for kk =(mm+1:numel(mocapstructcellarray))
    for jj = fieldsubset%:numel(comparefields)
output_struct{mm,kk,jj} = compare_features(mocapstructcellarray{mm},mocapstructcellarray{mm}.modular_cluster_properties,...
    mocapstructcellarray{kk},mocapstructcellarray{kk}.modular_cluster_properties,comparison_type,...
    mocapstructcellarray{mm}.annotated_pose_struct.(comparefields{jj}),...
    mocapstructcellarray{kk}.annotated_pose_struct.(comparefields{jj}),2,0);
kfoldlosses(mm,kk,jj) = output_struct{mm,kk,jj}.kfoldloss;
kfoldlosses(kk,mm,jj) = kfoldlosses(mm,kk,jj) ;

end
end
end


for mm =1:numel(mocapstructcellarray)
for kk =(mm+1:numel(mocapstructcellarray))
    for jj = fieldsubset%:numel(c
kfoldlosses(kk,mm,jj) = kfoldlosses(mm,kk,jj) 
    end
end
end
%% compare pre to the rest
figure(455)
imagesc(squeeze(kfoldlosses(1,2:end,fieldsubset ))')
set(gca,'XTick',1:numel(mocapstructcellarray)-1,'XTickLabels',names(2:end))
set(gca,'YTick',1:numel(fieldsubset ),'YTickLabels',annotation_fields(fieldsubset))
colorbar;
disp('all done \n')
caxis([0 0.25]);
title('comparison of pre lesion to full on the trunk')


figure(456)
imagesc(squeeze(kfoldlosses(2,:,fieldsubset ))')
set(gca,'XTick',1:numel(mocapstructcellarray),'XTickLabels',names(1:end))
set(gca,'YTick',1:numel(fieldsubset ),'YTickLabels',annotation_fields(fieldsubset))
colorbar;
disp('all done \n')
title('comparison of pre lesion to full on the trunk')


for mm=fieldsubset 
    figure(460+mm)
    imagesc(squeeze(kfoldlosses(:,:,mm )))
set(gca,'XTick',1:numel(mocapstructcellarray),'XTickLabels',names(1:end))
set(gca,'YTick',1:numel(mocapstructcellarray),'YTickLabels',names(1:end))
title(annotation_fields{mm})
c =colorbar;
caxis([0 0.25]);
c.Label.String = 'kFold loss';
end

fieldsubset = [1,3:6];
markersubset = [1,2,7,9];
%markersubset = [3:6];

kfoldlosses_markers = zeros(numel(annotation_fields),numel(subcluster_names));

mm =1;
kk =2;

    for jj = fieldsubset%:numel(comparefields)
        for ll = markersubset
output_struct_temp = compare_features(mocapstructcellarray{mm},mocapstructcellarray{mm}.modular_cluster_properties,...
    mocapstructcellarray{kk},mocapstructcellarray{kk}.modular_cluster_properties,comparison_type,...
    mocapstructcellarray{mm}.annotated_pose_struct.(comparefields{jj}),...
    mocapstructcellarray{kk}.annotated_pose_struct.(comparefields{jj}),ll,0);
kfoldlosses_markers(jj,ll) = output_struct_temp.kfoldloss;

end
    end

%% compare pre to the rest
figure(555)
imagesc(squeeze(kfoldlosses_markers(fieldsubset,markersubset)))
set(gca,'XTick',1:numel(markersubset),'XTickLabels',subcluster_names(markersubset))
set(gca,'YTick',1:numel(fieldsubset ),'YTickLabels',annotation_fields(fieldsubset))
colorbar;
disp('all done \n')
title('comparison of pre lesion to full on the trunk')

    
    


end