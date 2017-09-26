%% create all masks for a set of directories
function createmasks(camera_base_agg,numcameras)



for klk = 1:numel(camera_base_agg)
    if (numel(camera_base_agg{klk}))
    %% turn this into a separate modular file
    for cameratype_ind = 1:numcameras
        suffix = {'U','L','R'};
        cameradirectory = strcat(camera_base_agg{klk},filesep,'Camera',suffix{cameratype_ind});
        imageseries_folder = cameradirectory;
        
        mask_filename = strcat(cameradirectory,filesep,'mask_vals.mat');
        if (~exist(mask_filename,'file'))
            
            %% conver the mkv files to mp4 files
            fprintf('converting mkv files \n')
            read_mkv_file(imageseries_folder,cameradirectory );
            
            %imds = imageDatastore(moviefiles(3).name,'ReadFcn',@mpfour_reader);
            
            moviefiles = dir(strcat(cameradirectory,filesep,'*.mkv'));
            videofeatures = [];
            filenames = zeros(1,numel(moviefiles));
            for ll =1:numel(moviefiles)
                filenames(ll) = str2num(strrep(moviefiles(ll).name,'.mkv',''));
            end
            [vals,ind] = sort( filenames);
            fprintf('number of movie files %f \n',numel(moviefiles));
            
            ll=1
            movie_name_here = strcat(cameradirectory,filesep,num2str(vals(ll)),'.mkv');
            fprintf('For videofile %s loading in images \n',movie_name_here)
            tic
            videofile = [];
            %  imagefile = strcat(cameradirectory,filesep,movie_name_here);
            temp = mpfour_reader(movie_name_here,1);
            videofile = cat(2,videofile,temp);
            
            videomatrix = single(cat(4,videofile(:).cdata));
            
            mean_videomatrix = mean(videomatrix,4);
            
            videomatrixnormalized = bsxfun(@minus,videomatrix,mean_videomatrix);
            
            figure(282)
            imagesc(squeeze(mean(videomatrix(:,:,:,1),3)))
            mask_vals = imrect;
            mask_position = mask_vals.getPosition;
            
            videomatrixnormalized_roi = videomatrixnormalized(mask_position(2):(mask_position(2)+mask_position(4)),...
                mask_position(1):(mask_position(1)+mask_position(3)),:,:);
            
            figure(283)
            imagesc(squeeze(mean(videomatrixnormalized_roi(:,:,:,1),3)))
            
            save( mask_filename,'mask_position')
            
            
        end
    end
    end
end
end