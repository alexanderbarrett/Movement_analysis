function makeXaxislabel(NumTicks,ticklabels)

%xlim([1 NumTicks])
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
xlabel('day')
%legend('d1','d2','both')
set(gca,'XTickLabel', ticklabels)
end