function makeaxislabel(NumTicks,ticklabels)

ylim([1 NumTicks])
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
xlabel('day')
%legend('d1','d2','both')
set(gca,'YTickLabel', ticklabels)
xlabel('time after last lick')
end