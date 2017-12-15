function plot_intermarker_distribution(mocapstruct,time_subset,colorin)


 %% break this up into its own script
 
 
 
 %% knee up/down diff
 figure(500)
 
 subplot(1,4,1)
 hold on
 [n,x] = hist(( mocapstruct.markers_aligned_preproc.KneeL(time_subset,1)-mocapstruct.markers_aligned_preproc.KneeR(time_subset,1)),-100:1:100);
b1 =bar(x,n./sum(n),'EdgeColor','none','FaceColor',colorin)
%b1 = bar(x2,n2./sum(n2),'r','EdgeColor','none')
b1.FaceAlpha = 0.5;
title('Knee L-R Y(along axis)')

 subplot(1,4,2)
 hold on
 [n,x] = hist(( mocapstruct.markers_aligned_preproc.HeadF(time_subset,2)-mocapstruct.markers_aligned_preproc.HeadF(time_subset,2)),-100:1:100);
b1 =bar(x,n./sum(n),'EdgeColor','none','FaceColor',colorin)
%b1 = bar(x,n./sum(n2),'r','EdgeColor','none')
b1.FaceAlpha = 0.5;
title('HeadF-SpineF X (transverse axis)')

subplot(1,4,3)
 hold on
 [n,x] = hist(( mocapstruct.markers_aligned_preproc.SpineF(time_subset,3)-mocapstruct.markers_aligned_preproc.SpineM(time_subset,3)),-100:1:100);
b1 =bar(x,n./sum(n),'EdgeColor','none','FaceColor',colorin)
%b1 = bar(x2,n2./sum(n2),'r','EdgeColor','none')
b1.FaceAlpha = 0.5;
title('SpineF-SpineM Z')

subplot(1,4,4)
 hold on
 [n,x] = hist(( mocapstruct.markers_aligned_preproc.SpineM(time_subset,2)-mocapstruct.markers_aligned_preproc.SpineL(time_subset,2)),-100:1:100);
b1 =bar(x,n./sum(n),'EdgeColor','none','FaceColor',colorin)
%b1 = bar(x2,n2./sum(n2),'r','EdgeColor','none')
b1.FaceAlpha = 0.5;
title('SpineM-SpineL (along axis)')
 


%others: shouylder L and R

end