function [av_vel,vel_comps,av_std,std_comps] = get_markersubgroup_velocity(markermocap,submarkers,params)

fieldnames_here = fieldnames(markermocap);

%velocity
av_vel = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
vel_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);

%standard deviation
av_std = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
std_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);
for jj = 1:3
    vel_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
        std_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);

    for mm = submarkers
        [veloc,stand] = get_vector_velocity(markermocap.(fieldnames_here{(mm)})(:,jj),params);
vel_temp = vel_temp +veloc./numel(submarkers);
std_temp = std_temp +stand./numel(submarkers);

    end
    vel_comps(:,jj) = vel_temp;
    av_vel = av_vel + vel_temp.^2;
    
      std_comps(:,jj) = std_temp;
    av_std = av_std + std_temp.^2;
end
av_std = sqrt(av_std);
av_vel = sqrt(av_vel);

end
