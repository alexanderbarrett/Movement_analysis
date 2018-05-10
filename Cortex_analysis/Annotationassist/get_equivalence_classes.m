function [predicted_labels_equivalent,new_behavioral_labels] = get_equivalence_classes(predicted_labels_num,behavior_names)

    %define equivalence classes
    equivalences = {{'prone_still','StillProne'},...
    {'RearStill','rear_still'},...
    {'AdjustSteps','HeadRaise','HeadLower','HeadRotate','ArmAdjust'},...
    {'Anogenitalgroom','AnogenitalGroom'}...
    {'Swap', 'FspineSwap'} };

    renames = {'QuietRest','QuietRear','PosturalAdjustment','AnogenitalGroom','Swap'};
    % get numerical classes
    equivalence_ids = cell(1,numel(equivalences));
    for kk = 1:numel(equivalences)
        equivalence_ids{kk} = zeros(0,0);
        for ll = 1:numel(equivalences{kk})
            if find(strcmp(behavior_names,equivalences{kk}{ll}))
                equivalence_ids{kk} = cat(1, equivalence_ids{kk},find(strcmp(behavior_names,equivalences{kk}{ll})));
            else
                fprintf('string not found %f %f \n',kk,ll)
            end
        end
    end

%     %loop and set equal to first number/name
%     new_behavioral_labels = behavior_names;
%     predicted_labels_equivalent=predicted_labels_num;
%     for kk = 1:numel(equivalences)
%         predicted_labels_equivalent(ismember(predicted_labels_num,equivalence_ids{kk})) = equivalence_ids{kk}(1);
%         for nn = 1:numel(equivalence_ids{kk})
%             new_behavioral_labels{equivalence_ids{kk}(nn)} = renames{kk};
%         end
%     end

    %loop and set equal to first number/name
    new_behavioral_labels = behavior_names;
    predicted_labels_equivalent=predicted_labels_num;
    for kk = 1:numel(equivalences)
        predicted_labels_equivalent(ismember(predicted_labels_num,equivalence_ids{kk})) = equivalence_ids{kk}(1);
        
            new_behavioral_labels{equivalence_ids{kk}(1)} = renames{kk};
        
    end

end