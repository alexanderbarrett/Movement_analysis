function plot_score_histogram(valuesin,colorin,plotno)

axissize = size(valuesin);
[minval,minind] = min(axissize);
valuesin = permute(valuesin,[minind setxor([1,2],minind)]);
figure(1000+plotno)
axis_hist = -200:1:200;

for ind_pick = 1:minval
    subplot_tight(4,4,ind_pick)
 %   ind_pick
    
    [n,x] = hist( squeeze(valuesin(ind_pick,:)),axis_hist);
    n = n./sum(n);
  %  bar(x,log10(n),'b','EdgeColor','none');
       semilogy(x,(n),'Color',colorin,'Linewidth',2);
     hold on
end
%print('-dpng',strcat(plotdir_here,'posture_eigenvalue_hists.png'))
end
