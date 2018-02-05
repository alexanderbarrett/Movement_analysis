function [transout,statefreq] = compute_transition_matrix(annotation_vector,fieldnames)
maxval =max(annotation_vector); %include null/rest as 0

transprob = zeros(maxval+1,maxval+1);


transition_state = find(strcmp(fieldnames,'Transition'));
% remove states marked 'transition'
if transition_state
annotation_vector((annotation_vector == (transition_state-1))) = [];
end

statefreq = arrayfun(@(x) numel(find(annotation_vector == x)),0:maxval);
nz_states = find(statefreq>0);


for mm = 1:numel(annotation_vector)-1
    transprob(annotation_vector(mm)+1,annotation_vector(mm+1)+1) = transprob(annotation_vector(mm)+1,annotation_vector(mm+1)+1)+...
        1./statefreq(annotation_vector(mm)+1);
end
transprob_offdiag = transprob-diag(diag(transprob));
%normalize by the size of the diagonal
exit_prob = bsxfun(@rdivide,transprob_offdiag,(1-diag(transprob)));

figure(77)
imagesc(exit_prob (nz_states,nz_states));
set(gca,'XTick',1:numel(nz_states),'XTickLabels',fieldnames(nz_states))
set(gca,'YTick',1:numel(nz_states),'YTickLabels',fieldnames(nz_states))
xtickangle(90)
colorbar

figure(79)
bar(diag(transprob));
set(gca,'XTick',1:numel(nz_states),'XTickLabels',fieldnames(nz_states))
xtickangle(90)
title('self transition prob')
colorbar

figure(78)
bar(statefreq(nz_states)./numel(annotation_vector));
set(gca,'XTick',1:numel(nz_states),'XTickLabels',fieldnames(nz_states))
xtickangle(90)
title('state frequency')

%% get length dist

    cluster_size_dist = cell(1,max(annotation_vector));
    for ll = 1:max(annotation_vector)+1
        inst_label = zeros(1,numel(annotation_vector));
        inst_label(find(annotation_vector == ll+1)) = 1;
        pixellist = bwconncomp(inst_label);
        cluster_size_dist{ll} = cellfun(@numel,pixellist.PixelIdxList);        
    end
    
    fractionshort = cellfun(@(x) numel(x(x<=2)),cluster_size_dist)./cellfun(@numel,cluster_size_dist);

figure(80)
bar(fractionshort(nz_states));
set(gca,'XTick',1:numel(nz_states),'XTickLabels',fieldnames(nz_states))
xtickangle(90)
title('number of states<2 frames in length')
    
end