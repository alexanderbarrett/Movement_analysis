function av_vel = get_markersubgroup_velocity(markermocap,submarkers,params)

fieldnames_here = fieldnames(markermocap);

av_vel = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
for jj = 1:3
    vel_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
    for mm = submarkers
vel_temp = vel_temp +get_vector_velocity(markermocap.(fieldnames_here{(mm)})(:,jj),params)./numel(submarkers);
    end
    av_vel = av_vel + vel_temp.^2;
end
av_vel = sqrt(av_vel);

end
