function cellTraceInds = findExampleTracesIndividually(cellTraces, nTraces, cellColors, cellMapImage, centroids, cellTraceInds, varargin)

%%% Use this to get merged image and colors:
%%% [mergedImage, cellColors] = mergeImages(cellImages2);


nCellsTotal=size(cellTraces,1);
timeFrame=801:1200;

cellTraces=cellTraces(:,timeFrame);
maxVals=max(cellTraces,[],2);
[~,maxInds]=sort(maxVals,'descend');
cellTraces=cellTraces(maxInds,:);
maxVal=prctile(maxVals,90);

if isempty(cellTraceInds)
    cellTraceInds=zeros(nTraces,1);
    cellInd=1;
    saveTraceInd=1;
    while cellInd<nCellsTotal && saveTraceInd<nTraces+1

        h=figure(87);
        set(h, 'CurrentCharacter', 'k')
        plot(cellTraces(cellInd,:), 'Linewidth', 2, 'Color', cellColors(cellInd,:));
        ylim([0 max(maxVal,max(cellTraces(cellInd,:)))]); title([num2str(cellInd) ' of ' num2str(nCellsTotal)])
        waitforbuttonpress;
        reply=get(h, 'CurrentCharacter');

        switch reply
            case 'y'
                cellTraceInds(saveTraceInd)=cellInd;
                saveTraceInd=saveTraceInd+1;
            case 'q'
                cellInd=nCellsTotal;
            case 'b'
                cellInd=cellInd-2;
        end


        cellInd=cellInd+1;
    end
end    

figure;
for traceInd=1:length(cellTraceInds)
    cellInd=cellTraceInds(traceInd);
    if cellInd>0
        plot(cellTraces(cellInd,:)-(traceInd-1)*0.8*maxVal, 'Linewidth', 2, 'Color', cellColors(cellInd,:))
        hold on;
    end
end
set(gca, 'XTick', [], 'YTick', []); box off;

if ~isempty(cellMapImage) && ~isempty(centroids)
    figure; imagesc(cellMapImage); hold on; plot(centroids(cellTraceInds,1), centroids(cellTraceInds,2), 'w.', 'Markersize', 12)
    strValues = strtrim(cellstr(num2str((1:nTraces)','%d ')));
    text(centroids(cellTraceInds,1)+2, centroids(cellTraceInds,2)+2, strValues, 'Color', [1 1 1], 'Fontsize', 14)
    set(gca, 'XTick', [], 'YTick', [])
end