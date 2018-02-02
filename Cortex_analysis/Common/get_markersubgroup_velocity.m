function [av_vel,vel_comps,av_std,std_comps,av_accel,accel_comps,av_std_accel,std_comps_accel] = get_markersubgroup_velocity(markermocap,submarkers,params)

fieldnames_here = fieldnames(markermocap);

%velocity
av_vel = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
vel_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);

%standard deviation
av_std = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
std_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);

%accel and accel std
av_accel = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
accel_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);

av_std_accel = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
std_comps_accel = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);


for jj = 1:size(markermocap.(fieldnames_here{(1)}),2)
    vel_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
        std_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
         accel_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);
        accel_std_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),1);

    for mm = submarkers
        [veloc,stand,accelout,accelstd] = get_vector_velocity(markermocap.(fieldnames_here{(mm)})(:,jj),params);
vel_temp = vel_temp +veloc./numel(submarkers);
std_temp = std_temp +stand./numel(submarkers);

accel_temp = accel_temp +accelout./numel(submarkers);
av_std_accel = accel_std_temp +accelstd./numel(submarkers);

    end
    vel_comps(:,jj) = vel_temp;
    av_vel = av_vel + vel_temp.^2;
    
      std_comps(:,jj) = std_temp;
    av_std = av_std + std_temp.^2;
    
    %% accel
 accel_comps(:,jj) = accel_temp;
    av_accel = av_accel + accel_temp.^2;
    
      std_comps_accel(:,jj) = av_std_accel;
    av_std_accel = av_std_accel + av_std_accel.^2;
end
av_std = sqrt(av_std);
av_vel = sqrt(av_vel);

%%accel
av_std_accel = sqrt(av_std_accel);
av_accel = sqrt(av_accel);
end
