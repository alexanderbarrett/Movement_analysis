function clustmembers = entplot(constitcells, totent, entthresh, fntitle)
% Schnitzer & Meister, Fig. 2B

subplot(2,1,1)
% nICs = length(unique(constitcells(:,[1,2])))-1;
nICs = size(constitcells,1);
nICsorig = nnz(constitcells(:,2)==0);
colord = 'rb';
for j=1:nICs
    plot([0, totent(j)], j*[1,1], '-', 'Color', colord(1+(j<=nICsorig)),'LineWidth',2)
    hold on
    if j>nICsorig
        plot([1,1,1]*totent(j), [j,constitcells(j,1:2)], 'k-','LineWidth',2)
        plot([1,1]*totent(j), [constitcells(j,1:2)], 'ko','MarkerFaceColor','k','LineWidth',2)
    end
end
hold off
set(gca,'YDir','reverse','XDir','reverse')
set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','i','LineWidth',2,'TickLength',[0.01,0.02])
xlabel('Entropy of compressed spike train')
ylabel('Cell cluster number')
xlim([totent(end)*.9999, totent(1)*1.0001])
ylim([0,nICs+1])
set(gcf,'Color','w','PaperPositionMode','auto')
title(['\DeltaH threshold = ',num2str(entthresh,3)])

% Determine the original cells corresponding to final clusters
subplot(2,1,2)
for j=nICs:-1:1
    clustmembers{j} = constitcells(j,1:2);
    cellgrps = clustmembers{j}(clustmembers{j}>nICsorig);
    while ~isempty(cellgrps)
        cnew = clustmembers{j}(clustmembers{j}<=nICsorig);
        for k = cellgrps
            cnew = [cnew,constitcells(k,1:2)];
        end
        clustmembers{j} = cnew(cnew>0);
        cellgrps = clustmembers{j}(clustmembers{j}>nICsorig);
    end
    clustmembers{j} = clustmembers{j}(clustmembers{j}>0);
    fprintf('Cluster %3.0f contains %3.0f cells.\n', j, length(clustmembers{j}))
end
colord = 'krgbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb';
symbord = 'osvdxphhhhhhhhhhhhhhhhhhhhhhhhhh';
hlabel = [];
for j=nICsorig+1:nICs
    hlabel(length(clustmembers{j})-1) = plot(j*ones(1,length(clustmembers{j})) - nICsorig, clustmembers{j}, ...
        [symbord(length(clustmembers{j})), '-'], 'MarkerSize', 10, ...
        'Color', colord(length(clustmembers{j})), 'MarkerFaceColor', colord(length(clustmembers{j})), 'LineWidth', 2);
    hold on
end
hold off
set(gca,'FontSize',22,'FontWeight','bold','FontName','Helvetica','FontAngle','i','LineWidth',2,'TickLength',[0.01,0.02])
xlabel('Cluster number')
ylabel('Original cells participating')
ylim([0,nICsorig+1])
xlim([0,nICs-nICsorig+1])
legend(hlabel', {'Cell pairs','Triples','4 cells','5 cells','>=6 cells'}, ...
    'Location','EastOutside')
set(gcf,'Color','w','PaperPositionMode','auto')
if nargin>=5
    title(fntitle)
end

% subplot(2,1,2)
% dendmat = constitcells(nICsorig+1:end,:);
% dendmat(:,3) = cumsum(constitcells(nICsorig+1:end,3));
% % dendmat = constitcells;
% % dendmat(:,3) = cumsum(constitcells(:,3));
% dendmat = flipud(dendmat)
% dendrogram(dendmat)
