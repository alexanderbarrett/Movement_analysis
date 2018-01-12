function plot_behavior_examples(mocapstruct,fulloutput,allbouts)





params.fps = 300;
markerout = hipass_clip_fragments(mocapstruct.markers_aligned_preproc ,fulloutput,params);

xlimhere = 5000;

markerstoplot = {'HeadF','Offset1','ElbowL','ElbowR','KneeR','KneeL'};
numplots = numel(markerstoplot);
boutstoplot = min(5,numel(allbouts));

%% plot PSD and time trace examples
for kk = 1:numel(markerstoplot)
    
    figure(77)
subplot(1,numplots,kk)
for ll = 1:boutstoplot
    indshere =   arrayfun(@(x) find(fulloutput == allbouts{ll}(x)),1:numel(allbouts{ll}));

plot(markerout.(markerstoplot{kk})(indshere,3));
hold on
end
title(markerstoplot{kk})



  
    figure(88)
subplot(1,numplots,kk)

for mm = 1:3
 [Pxx2,fxx] = pwelch(markerout.(markerstoplot{kk})(:,mm),mocapstruct.fps*2,floor(0.5*mocapstruct.fps),mocapstruct.fps*2,mocapstruct.fps,'onesided');

 if (mm == 1)
     P_sum2 = Pxx2;
 else
    P_sum2 = P_sum2 + Pxx2; 
  end
 end
% 
% 
 plot(fxx,log10(P_sum2),'b')

 xlim([1 50])
 ylabel('SP ')
 xlabel('Freq')
 title(markerstoplot{kk}) 
end


%% plot animation subset
h = figure(370)
marker_inds = 1:20;
agg_features = [];
bout_use = 2;
for mm =marker_inds
agg_features = cat(3,agg_features,mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mm})(allbouts{bout_use} ,:));
end
pictoral_animation_subset(permute(agg_features,[1 3 2]),mocapstruct,(1:10:numel(allbouts{bout_use})),175,marker_inds,h)
xlim([0 1500])


end