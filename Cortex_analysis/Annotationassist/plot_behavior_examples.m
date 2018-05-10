function plot_behavior_examples(mocapstruct,fulloutput,allbouts)

params.fps = 300;
markerout = hipass_clip_fragments(mocapstruct.markers_aligned_preproc ,fulloutput,params);

xlimhere = 5000;

markerstoplot = {'HeadF','Offset1','ElbowL','ElbowR','KneeR','KneeL'};
numplots = numel(markerstoplot);
boutstoplot = min(1,numel(allbouts));
[~,bouttoplot] = max(cellfun(@numel,allbouts));

%% plot PSD and time trace examples
legendnames = cell(1,0);
for kk = 1:numel(markerstoplot)
    
    
    %% individual subplots (for comparing bouts)
    figure(77)
subplot(1,numplots,kk)
for ll = bouttoplot
    indshere =   arrayfun(@(x) find(fulloutput == allbouts{ll}(x)),1:numel(allbouts{ll}));

plot(markerout.(markerstoplot{kk})(indshere,3));
hold on
end
title(markerstoplot{kk})

%% all on one plot
figure(78)
colorlist = jet(numel(markerstoplot));
 indshere2 = [];
 breakpoints = [];
for ll = [bouttoplot]
    indshere2 =   cat(2,indshere2,arrayfun(@(x) find(fulloutput == allbouts{ll}(x)),1:numel(allbouts{ll})));
    breakpoints = cat(2,breakpoints,numel(indshere2));
end
plot(1000*(0:1./params.fps:((numel(indshere2)-1)./params.fps)),30*(kk-1)+markerout.(markerstoplot{kk})(indshere2,3),'Linewidth',2,'Color',colorlist(kk,:));
hold on

legendnames{kk} = markerstoplot{kk};

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
figure(78)
for zz = 1:numel(breakpoints)
       plot(1000*repmat(breakpoints(zz),1,12)./params.fps-1./params.fps,(-20*2):20:9*20,'k','linewidth',2) 
end
legend(legendnames,'Location','EastOutside')
xlabel('Time (ms)')
ylabel('Relative Marker z-position')
box off

%% plot animation subset
h = figure(370)
marker_inds = 1:20;
agg_features = [];
bout_use = 1;
for mm =marker_inds
agg_features = cat(3,agg_features,mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mm})(allbouts{bout_use} ,:));
end
pictoral_animation_subset(permute(agg_features,[1 3 2]),mocapstruct,(1:10:numel(allbouts{bout_use})),175,marker_inds,h)
xlim([0 1500])


end