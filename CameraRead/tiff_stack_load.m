function output_movie = tiff_stack_load(directory,camno,idstring)
search_directory = strcat(directory,'*CAMERA',num2str(camno),'*',idstring,'*');
tiff_names = dir(search_directory);

if (numel(tiff_names) == 0)
fprintf(strcat('No tiff files by that name found in  ',search_directory,' \n')   )
end


[image_1,map] = imread(strcat(directory,tiff_names(1).name));
if (numel(map) == 0)
    [~, map] = gray2ind(image_1, 256);
end

%test=struct();

for k = 1:numel(tiff_names)
[image_1] = imread(strcat(directory,tiff_names(k).name));
im_indexed= gray2ind(image_1, 256);
output_movie(k) = im2frame(im_indexed,map);
end
 movie(output_movie);


%deal with colormap and formating


%immat = zeros(
%output_tiff = 
%read

end