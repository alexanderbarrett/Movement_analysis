%% master script for viewing movies
%Jesse Marshall 20160726

file_directory = 'E:/Bence/Data/MOCAP/';
file_directory_day = strcat(file_directory,'20160725/');
tiff_cam1 = tiff_stack_load(file_directory_day,1,'move5');