function [D] = find_weighted_distance(data)
%finds the weighted distances (D) between all rows in 'data'

    N = length(data(:,1));
    
    %want a(a-b)
    D= repmat(nansum(data.^2,2),1,N)-data*data';
   % D= pdist(data);
   % logData = log(data);
    %logData(isinf(logData) | isnan(logData)) = 0;
    
   % entropies = -sum(data.*logData,2);
    
   % D = - data * logData';
    
   % D = bsxfun(@minus,D,entropies);
    
   % D = D ./ log(2);
   % D(1:(N+1):end) = 0;