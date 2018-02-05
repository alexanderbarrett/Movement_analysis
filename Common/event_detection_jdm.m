%taken from Laurie's code
%jesse marshall oct 2012

function [finalpeaks,offsetpeaks,peakrise,thresh] = event_detection_jdm(tracesfile,imagesfile,params,raster_file)

close all;

load(tracesfile);
load(imagesfile)

%% for Session2 mice (980 and 400 written with the new code)
%  if (strcmp(class(SpikeTraceData),'double'))
%         temp = SpikeTraceData;
%        SpikeTraceData = struct([]);
%         for jj = 1:size(temp,2)
%             SpikeTraceData(jj).Trace = temp(:,jj);
%         end        
%  end
 
 if (exist('IcaTraces','var'))
         SpikeTraceData = struct([]);
      for jj = 1:size(IcaTraces,1)
    SpikeTraceData(jj).Trace = squeeze(IcaTraces(jj,:,:))';
      end
 end

 if (exist('IcaFilters','var'))
         SpikeImageData = struct([]);
      for jj = 1:size(IcaFilters,1)
    SpikeImageData(jj).Image = squeeze(IcaFilters(jj,:,:))';
      end
 end
 
%% Declare vars
doOrdFilt = 0;%get(handles.DoOrdFilt, 'Value');
ordFiltOrder = 1; %str2double(get(handles.OrdFiltOrder, 'String'));
ordFiltDomain = 1; %str2double(get(handles.OrdFiltDomain, 'String'));
makeThePlot = params.makeThePlot; %get(handles.MakeThePlot, 'Value');
doMovAvg = params.doMovAvg;%get(handles.DoMovingAverage, 'Value');
movAvgFiltSize = params.NumFramesLocalAverage; % str2double(get(handles.NumFramesLocalAverage, 'String'));
numStdsForThresh = params.numStdsForThresh;% str2double(get(handles.StdDevs, 'String'));
movAvgReqSize = params.movAvgReqSize;% str2double(get(handles.MovAvgReqSize, 'String'));
reportMidpoint = 1; %get(handles.ReportMidpoint, 'Value');
minTimeBtEvents=  1; %what units? str2double(get(handles.MinTimeBtEvents,'String'));


%for neighbor routine - hacked the object a bit
handles.OverlapRadius = 10; %microns?

%% Pull Out traces
tracesToProcess=1:length(SpikeTraceData);
waitbarInterval=round(length(tracesToProcess)/10);
traceLength=length(SpikeTraceData(tracesToProcess(1)).Trace);
celltraces=zeros(length(tracesToProcess), traceLength);


for ind=1:length(tracesToProcess)
    trInd=tracesToProcess(ind);
    celltraces(ind, :)=SpikeTraceData(trInd).Trace;
    
end

imagesOfICs=1:length(SpikeImageData);
icSize=size(SpikeImageData(imagesOfICs(1)).Image);
icmat=zeros([icSize, length(imagesOfICs)]);

for ind=1:length(imagesOfICs)
    icInd=imagesOfICs(ind);
    icmat(:,:,ind)=SpikeImageData(icInd).Image;
end



%use um space between criteria

neighborsCellFull = identifyNeighborsAuto(icmat, handles);
fprintf('Done assigning neighbor IDs ...')



%% Set up the analysis and calculate the SD of the full traces set
numCells = size(neighborsCellFull,1);
cellstdAll = zeros(numCells,1);

numCells


% Get SD by fitting gaussian to the histogrammed data (because of heavy pos
% tail from bursting)
for cellnum = 1 : numCells
    xdata = linspace(-0.3,0.3,1000)';         %%%% might want to take a look at these values
    if doOrdFilt
        ydata = hist(celltraces(cellnum,:) - ...
            ordfilt2(celltraces(cellnum,:)',ordFiltOrder,ones(ordFiltDomain,1),'symmetric')',xdata)';
    else
        ydata=hist(celltraces(cellnum,:)-mean(celltraces(cellnum,:)), xdata)';
    end
    
    options = fitoptions('method','nonlinearleastsquares',...
        'startpoint',[100 0 0.01]);
    try
        fitres = fit(xdata,ydata,'gauss1',options);
        cellstdAll(cellnum,1) = (fitres.c1/sqrt(2));
    catch
        disp('fit failed...')
        cellstdAll(cellnum,1) = std(celltraces(cellnum,:));
    end
end
% clear celltraces
clear options fitres xdata ydata cellnum plotcount plottingOn

fprintf('Done calculating trace SD. \n')



%% Run the peak finding
%superiorPeaks = zeros(0,5);
origpeaks = cell(numCells,1);
finalpeaks = cell(numCells,1);
offsetpeaks = cell(numCells,1);
peakdecay = cell(numCells,1);
peakrise = cell(numCells,1);
thresh = zeros(numCells,1);
peak_amplitudes = cell(numCells);
peak_times = cell(numCells);


% get the traces shifted by the filtered trace
if doOrdFilt
    inputtraces = celltraces'  - ordfilt2(celltraces', ordFiltOrder,ones(ordFiltDomain,1),'symmetric');
    inputtraces = inputtraces';
else
    inputtraces=celltraces;
end


for c=1:size(celltraces,1)
    offsetpeaks{c} = zeros(1,0);
    % set the threshold to a multiple of the SD of the trace for that cell
    thresh(c) = numStdsForThresh*cellstdAll(c);
    
    if doMovAvg
        inputsignal = filtfilt(ones(1,movAvgFiltSize)/movAvgFiltSize,1,inputtraces(c,:));
        %         inputsignal = 1/2*(medfilt1(inputtraces(c,:),movAvgFiltSize)+...       %%% can remove or play with this avg/med combination
        %             filtfilt(ones(1,movAvgFiltSize)/movAvgFiltSize,1,inputtraces(c,:)));
    else
        inputsignal=inputtraces(c,:);
    end
    
    [~,testpeaks] = findpeaks(inputsignal,'minpeakheight',thresh(c));
    [~,testpeaks2] = findpeaks(inputsignal,'minpeakdistance',minTimeBtEvents);
    testpeaks = intersect(testpeaks,testpeaks2);
    clear testpeaks2
    testpeaks = intersect(testpeaks,find(filtfilt(ones(1,movAvgReqSize)/movAvgReqSize,1,...
        inputtraces(c,:))>thresh(c)));
    
    
    if isempty(testpeaks)
        origpeaks{c} = zeros(1,0);
        finalpeaks{c} = zeros(1,0);
        peakrise{c} = zeros(1,0);
        offsetpeaks{c} = zeros(1,0);
        continue
    end
    if isempty(neighborsCellFull{c,1})
        [vectowrite,vectoramplitudes] = findRecentTroughs(inputsignal,inputsignal,testpeaks);
        okpeaks = vectoramplitudes > thresh(c); % get the peaks with increase greater than thresh
        
        origpeaks{c} = testpeaks(okpeaks);
        finalpeaks{c} = testpeaks(okpeaks);
        if reportMidpoint
            offsetpeaks{c} = 1/2*(testpeaks(okpeaks) + vectowrite(okpeaks));
        else
            offsetpeaks{c} = testpeaks{okpeaks};
        end
        peakrise{c} = vectoramplitudes(okpeaks);
    else
        [vectowrite,vectoramplitudes] = findRecentTroughs(inputsignal,inputsignal,testpeaks);
        okpeaks = vectoramplitudes > thresh(c);
        mean_amplitudes = mean(okpeaks);
        % don't write into 'finalpeaks'
        origpeaks{c} = testpeaks(okpeaks);
        % %             finalpeaks{c,m} = testpeaks(okpeaks);
        if reportMidpoint
            offsetpeaks{c} = 1/2*(testpeaks(okpeaks) + vectowrite(okpeaks));
        else
            offsetpeaks{c} = testpeaks{okpeaks};
        end
        peakrise{c} = vectoramplitudes(okpeaks);
    end
end

for c=1:size(celltraces,1)
    testpeaks = origpeaks{c};
    
    if doMovAvg
        filteredtrace = filtfilt(ones(1,movAvgFiltSize)/movAvgFiltSize,1,inputtraces(c,:));
        %         filteredtrace = 1/2*(filtfilt(ones(1,movAvgFiltSize)/movAvgFiltSize,1,inputtraces(c,:))+...    %%% can remove or play with this avg/med combination
        %             medfilt1(inputtraces(c,:),movAvgFiltSize));
    else
        filteredtrace=inputtraces(c,:);
    end
    % skip step if empty or no neighbors
    if isempty(testpeaks)
        continue
    end
    if isempty(neighborsCellFull{c,1})
        continue
    end
    
    
    
    
    % check each peak to see if (1) neighbor spikes w/in 2 frames &&
    % (2) neighbor amplitude is larger
    othersgreater = zeros(length(testpeaks),1);
    
    
    
    %% For now, ignore neighbor filtration
    %
    %     testcount = 0;
    %     for p = testpeaks
    %         testcount = testcount + 1;
    %         neighborcellampl = inputtraces(neighborsCellFull{c,1},p) > inputtraces(c,p);
    %         neighborcount = 0;
    %         neighborcellpeak = [];
    %         for n = [neighborsCellFull{c,1}]'
    %             neighborcount = neighborcount + 1;
    %             neighborcellpeak(neighborcount) = ~isempty(intersect(origpeaks{n},p-2:1:p+2));
    %         end
    %
    %         rankingmatrix = [];
    %         if any(size(neighborcellpeak)~=size(neighborcellampl))
    %             neighborcellpeak=neighborcellpeak';
    %         end
    %         if any(neighborcellampl & neighborcellpeak)
    %             othersgreater(testcount) = 1;
    %             %%% ranking matrix keeps information about all the neighbors,
    %             %%% whether they have a peak, peak size, cell number, etc
    %             %rankingmatrix(:,1) = neighborcellpeak; % 1/0 if have peak
    %             %rankingmatrix(:,2) = inputtraces(neighborsCellFull{c,1},p); % max amplitude around peak
    %             %rankingmatrix(:,3) = neighborsCellFull{c,1}; % cell number
    %
    %             %rankingmatrix = sortrows(rankingmatrix,[-1 -2]);        %%% sort by whether have peak, then by peak size
    %             %rankingmatrix = rankingmatrix(1,:);                     %%% only take winning neighbor
    %             %rankingmatrix(1,4) = intersect(origpeaks{rankingmatrix(1,3)},p-2:1:p+2);
    %             %rankingmatrix(1,5) = 1; % add the trial number - %%% obsolete now
    %             %superiorPeaks = cat(1,superiorPeaks,rankingmatrix(1,:));
    %         end
    %
    %
    %         %%% not sure why this is here, since peak finding above had a
    %         %%% minimum distance between peaks
    % %         if any(filteredtrace(max(1,p-4):1:min(p+4,length(filteredtrace)))>filteredtrace(p))
    % %             othersgreater(testcount) = true;
    % %         end
    %     end
    %     clear neighborcellampl neighborcellpeak p testcount peakcount
    %
    %
    %
    %
    
    
    %%%% check the 'inputtraces' and take max within 3 frames back of
    %%%% it??
    
    %         othersgreater = max(inputtraces(neighborsCellFull{c,1},testpeaks),[],1) > ...
    %             inputtraces(c,testpeaks);
    %         testpeaks = testpeaks(~othersgreater);
    finalpeaks{c} = testpeaks(~othersgreater);
    offsetpeaks{c} = offsetpeaks{c}(~othersgreater);
    peakrise{c} = peakrise{c}(~othersgreater);
    
    % fix the dimensionality of the empty matrix
    if isempty(finalpeaks{c})
        finalpeaks{c} = zeros(1,0);
        offsetpeaks{c} = zeros(1,0);
        peakrise{c} = zeros(1,0);
    end
    
    
end

save(raster_file,'finalpeaks','offsetpeaks','peakrise');

if makeThePlot
    
    result = 100;
    while (result > 0)
        
        prompt = 'Enter IC number to Examine [0 for exit] ';
        result = input(prompt);
        
        if (result >0)
            c=result;
            figure(29)
%             clf;
%             subplot(2,1,1)
%             hold all
%             plot(inputtraces(neighborsCellFull{c,1},:)')
%             plot(inputtraces(c,:),'k','linewidth',2)
%             plot(testpeaks,inputtraces(c,testpeaks),'r.','markersize',20)
%             plot(finalpeaks{c},inputtraces(c,finalpeaks{c}),'b.','markersize',20)
%             plot(round(offsetpeaks{c}),inputtraces(c,round(offsetpeaks{c})),'c.','markersize',16)
%             plot(2*offsetpeaks{c} - finalpeaks{c},inputtraces(c,2*offsetpeaks{c} - finalpeaks{c}),'g.','markersize',16)
%             %plot(peakrise{c},inputtraces(c,peakrise{c}),'g.','markersize',16)
%             
%             title(sprintf('main cell number %d',c))
%             hold off
%             legend(num2str(cat(1,neighborsCellFull{c,1},c)))
          %  subplot(2,1,2)
            hold on
            if (movAvgFiltSize > 1)
            plot(filtfilt(ones(1,movAvgFiltSize)/movAvgFiltSize,1,inputtraces(c,:)),'k','linewidth',2)
            else
            plot(inputtraces(c,:),'c','linewidth',2)    
            end
            plot(testpeaks,inputtraces(c,testpeaks),'c.','markersize',20)
            plot(finalpeaks{c},inputtraces(c,finalpeaks{c}),'b.','markersize',20)
            plot(round(offsetpeaks{c}),inputtraces(c,round(offsetpeaks{c})),'r.','markersize',16)
            plot(2*offsetpeaks{c} - finalpeaks{c},inputtraces(c,2*offsetpeaks{c} - finalpeaks{c}),'g.','markersize',16)
            hold off
        end
    end
end








%not sure what handles are




