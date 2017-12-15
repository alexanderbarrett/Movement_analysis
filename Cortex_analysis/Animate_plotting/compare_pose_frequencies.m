function compare_pose_frequencies(varargin)

disp('Number of input arguments: '+ nargin)
    celldisp(varargin)

    labels = fieldnames(varargin{1});
    
    barvals = zeros(numel(labels),nargin-1);
for kk=1:nargin-1
    barvals(:,kk) = structfun(@numel,varargin{kk})./...
            nansum(structfun(@numel,varargin{kk}));

end

figure(344)
    bar(barvals)
    set(gca,'XTickLabels',labels)
    legend(varargin{end})
    title('overall frequency of pose')
    
    
figure(345)
    bar(bsxfun(@rdivide,barvals,barvals(:,1)))
    set(gca,'XTickLabels',labels)
    legend(varargin{end})
    title('normalized frequency of pose')
    
    
    
end