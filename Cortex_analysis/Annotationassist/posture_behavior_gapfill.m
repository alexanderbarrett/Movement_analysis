function [pose_struct_out,globalbehavior_struct_out] = posture_behavior_gapfill(pose_struct,globalbehavior_struct)
pose_struct_out = pose_struct;
globalbehavior_struct_out = globalbehavior_struct;

fields_copy_beh_to_posture = {'LBodyGroom','LArmGroom','LArmScratch','RBodyGroom','RArmGroom','RArmScratch','AnogenitalGroom'};
for mm = 1:numel(fields_copy_beh_to_posture)
    if isfield(globalbehavior_struct_out,fields_copy_beh_to_posture{mm})
pose_struct_out.(fields_copy_beh_to_posture{mm}) = globalbehavior_struct_out.(fields_copy_beh_to_posture{mm});
    end
end

    if isfield(globalbehavior_struct_out,'Still')

        if isfield(pose_struct_out,'ExtProne')
            globalbehavior_struct_out.prone_still = intersect(globalbehavior_struct_out.Still,cat(2,pose_struct_out.Prone,pose_struct_out.ExtProne));
        else 
            globalbehavior_struct_out.prone_still = intersect(globalbehavior_struct_out.Still,pose_struct_out.Prone);
        end
       
%%%%%%
% for JDM25 caff
% globalbehavior_struct_out.prone_still = intersect(globalbehavior_struct_out.Still,cat(2,pose_struct_out.Prone,pose_struct_out.ExtProne));

% for Vicon8 caff
% globalbehavior_struct_out.prone_still = intersect(globalbehavior_struct_out.Still,pose_struct_out.Prone);
%%%%%%%%%

globalbehavior_struct_out.rear_still = intersect(globalbehavior_struct_out.Still,cat(2,pose_struct_out.ExtRear,pose_struct_out.MidRear,pose_struct_out.LowRear));
    end
end