function [pVals MI] = Mutual_Information_Shuffled(x_trace,y_trace,numValues_Y,shuffle_flag,num_shufs)
%Shuffled MI, takes in a matrix of cells and traces, and a 'y_trace' that
%consists of multiple time traces to compare the x trace to, to see if the
%x trace is modulated differentially, as assessed by shuffling. 
%modified from Lacey Kitsch (maybe originally from Laurie Burns) by Jesse
%Marshall, 1/6/2014

s = RandStream.create('mt19937ar','seed',sum(100*clock));
h=waitbar(0, 'calculating shuffled MI...');


% get parameters from interface
yTrace=y_trace;
num_y = size(y_trace,2);

numValuesY=numValues_Y;
numShufs=num_shufs;
numXTraces=size(x_trace,2);
numXtraces=size(x_trace,2);

shufMIvals=zeros(numShufs, numXtraces,num_y);



% construct logical event matrix, X


numPoints=size(x_trace,1);
X=false(numXtraces, numPoints);
for i=1:numXtraces
    
    X(i,:)=logical(x_trace(:,i));
    
end



% construct discrete Y
%     oldYtrace=yTrace;
%     oldYtime=SpikeTraceData(yTraceInd).XVector;
%     xTime=timeTrace;
%     yTrace=zeros(size(timeTrace));
%     for timeInd=1:numPoints
%         [~, closestYtimeInd]=min(abs(oldYtime-xTime(timeInd)));
%         yTrace(timeInd)=oldYtrace(closestYtimeInd);
%     end

% divide y to have numValuesY discrete values, range 0 to
% numValuesY-1
Y = zeros(numPoints,num_y);
MI = zeros(numXtraces,num_y);
    logProbY=zeros(num_y,numValuesY);

for i=1:num_y
    Y(:,i) =round((squeeze(yTrace(:,i))-min(squeeze(yTrace(:,i))))*(numValuesY-1)/(max(squeeze(yTrace(:,i)))-min(squeeze(yTrace(:,i)))));
    
    
    %%%%% calculate actual MI
    % calculate probabilities
    logProbX0=log(sum(1-X,2)+1)-log(numPoints+2);
    logProbX1=log(sum(X,2)+1)-log(numPoints+2);
    logProbX0givenY=zeros(numXtraces, numValuesY);
    logProbX1givenY=zeros(numXtraces, numValuesY);
    logProbX0andY=zeros(numXtraces, numValuesY);
    logProbX1andY=zeros(numXtraces, numValuesY);
    
    for yVal=0:numValuesY-1
        logProbY(i,yVal+1)=log(sum(squeeze(Y(:,i))==yVal)+1)-log(numPoints+numValuesY);
        theseX=X(:,squeeze(Y(:,i))==yVal);
        logProbX0givenY(:,yVal+1)=log(sum(1-theseX,2)+1)-log(size(theseX,2)+2);
        logProbX1givenY(:, yVal+1)=log(sum(theseX,2)+1)-log(size(theseX,2)+2);
        logProbX0andY(:, yVal+1)=log(sum(1-theseX, 2)+1)-log(numPoints+2*numValuesY);
        logProbX1andY(:,yVal+1)=log(sum(theseX, 2)+1)-log(numPoints+2*numValuesY);
    end
    % calculate MI
    MI(:,i)=sum(exp(logProbX0andY).*(logProbX0andY-repmat(logProbX0,1,numValuesY)-repmat(squeeze(logProbY(i,:)),numXtraces,1)),2)+...
        sum(exp(logProbX1andY).*(logProbX1andY-repmat(logProbX1,1,numValuesY)-repmat(squeeze(logProbY(i,:)),numXtraces,1)),2);
end


if (shuffle_flag)
    %%%%%% calculate shuffled MI
    parfor shufInd=1:numShufs
        shufX=false(size(X));
        
        % shuffle cells
        for n=1:numXtraces
            shufX(n,:)=X(n,randperm(s,numPoints));
        end
        
        % smooth traces and binarize
        %         for n=1:size(shufX,1)
        %             shufX(n,:)=logical(conv(single(shufX(n,:)), gausswin(N,5/sigma), 'same'));
        %         end
        for i=1:num_y
            
            % calculate probabilities
            logProbX0givenY=zeros(numXtraces, numValuesY);
            logProbX1givenY=zeros(numXtraces, numValuesY);
            logProbX0andY=zeros(numXtraces, numValuesY);
            logProbX1andY=zeros(numXtraces, numValuesY);
            
            for yVal=0:numValuesY-1
                theseX=shufX(:,squeeze(Y(:,i)==yVal);
                logProbX0givenY(:,yVal+1)=log(sum(1-theseX,2)+1)-log(size(theseX,2)+2);
                logProbX1givenY(:, yVal+1)=log(sum(theseX,2)+1)-log(size(theseX,2)+2);
                logProbX0andY(:, yVal+1)=log(sum(1-theseX, 2)+1)-log(numPoints+2*numValuesY);
                logProbX1andY(:, yVal+1)=log(sum(theseX, 2)+1)-log(numPoints+2*numValuesY);
            end
            
            % calculate MI
            shufMIvals(shufInd, :,i)=sum(exp(logProbX0andY).*(logProbX0andY-repmat(logProbX0,1,numValuesY)-repmat(squeeze(logProbY(i,:)),numXtraces,1)),2)+...
                sum(exp(logProbX1andY).*(logProbX1andY-repmat(logProbX1,1,numValuesY)-repmat(squeeze(logProbY(i,:)),numXtraces,1)),2);
        end
    end
end


% get p-vals for all i -- the fraction of shuffled MI calues that are
% bigger than the measured -- ie if the randomized carries more
% information.
pVals=zeros(numXtraces,num_y);
for i=1:num_y
    for trInd=1:numXtraces
        pVals(trInd,i)=sum(squeeze(shufMIvals(:,trInd,i))>=squeeze(MI(trInd,i)))/numShufs;
    end
end
        
end


%



