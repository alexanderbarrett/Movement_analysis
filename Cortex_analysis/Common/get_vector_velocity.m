function [vectorout,vectorstdout,vectorout_accel,vectorstdout_accel]= get_vector_velocity(vectorin,params)


tracesmooth = medfilt1(vectorin,params.medfiltorder);
gfilter = fspecial('gaussian',[50 1], params.gaussorder);
tracesmoothed = convn(tracesmooth,gfilter,'same');

vectorout_vel = cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);

vectorout = movmean(vectorout_vel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints','fill');%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
vectorstdout = movstd(vectorout_vel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints','fill');%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);

a_diff_order = floor(params.difforder./2);
vectorout_vel_a = cat(1,zeros((a_diff_order)-2,1),tracesmoothed(a_diff_order:end,1) - tracesmoothed(1:end-a_diff_order+1),1);

vectorout_accel = cat(1,zeros((a_diff_order)-2,1),vectorout_vel_a(a_diff_order:end,1) - vectorout_vel_a(1:end-a_diff_order+1),1);

vectorout_accel = movmean(vectorout_accel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints','fill');%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
vectorstdout_accel = movstd(vectorout_accel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints','fill');%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);

end