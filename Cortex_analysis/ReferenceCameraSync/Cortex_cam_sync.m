
        cameradirectory = strcat(camera_base_agg{klk},filesep,'Camera',suffix{cameratype_ind});
        imageseries_folder = cameradirectory;
        
        
        
        % if (cameratype_ind == 1)
        %     cameradirectory = cameradirectory_U;
        %     imageseries_folder = imageseries_folder_U;
        % else
        %     cameradirectory = cameradirectory_L;
        %     imageseries_folder = imageseries_folder_L;
        % end
        %
        % if (~exist(cameradirectory))
        %     mkdir(cameradirectory)
        % end
        
        %% conver the mkv files to mp4 files
        fprintf('converting mkv files \n')
        camfolder = strcat(imageseries_folder,cameradirectory);
             
        camfolder = 'E:\Bence\Data\Motionanalysis_captures\Vicon8\20170820\Recording_day5_overnight\CameraL\636388255619839101\';

        read_mkv_file(camfolder,camfolder );
        lk = 1;
        videofeatures_directory{lk} = camfolder;
        %imds = imageDatastore(moviefiles(3).name,'ReadFcn',@mpfour_reader);
        
        %moviefiles = dir(strcat(cameradirectory,filesep,'*.mkv'));
        
        
        %.times file
                times_files = dir(strcat(videofeatures_directory{lk},filesep,'*.times'));

                 
        f = fopen(strcat(videofeatures_directory{lk},filesep,times_files.name));
        float1 = fread(f,[1,100000000],'uint64');
        frame_number = numel(float1);
        % 1000.*(frame_number./(float1(end)-float1(1)));
        frame_chunk = 5000;
        
                
                   %% resample at the variable framerate
        %true_video_framerate = fps*(original_length+450+video_missingframes*(video_fps./fps))./size(markers_preproc.(marker_names{1}),1);
        true_video_framerate = 59.82;
        % video_pc_traces_resampled{lk} = resample(double(video_pc_traces(1:num_pcs,:))',round(fps*100),round(true_video_framerate*100));
       
        %preallocate
        video_pc_traces_resampled{lk} = zeros( floor(frames_to_analyze*300./true_video_framerate*1.1),num_pcs);
        frames_total = 0;
        frame_current = 0;
        fprintf('resampling video files \n')
        for kk = 1:frame_chunk: frames_to_analyze
            chunk_start = kk;
            chunk_end = kk+frame_chunk-1;
            frame_chunk_here = frame_chunk;
            if (chunk_end> frames_to_analyze)
                chunk_end =  frames_to_analyze;
                frame_chunk_here = (chunk_end-chunk_start+1);
            end
            
            resample_rate = 1000.*(frame_chunk_here./(float1(chunk_end)-float1(chunk_start)));
            fprintf('chunk %f resample rate %f \n',kk,resample_rate);
            mult_factor = 200;
            resampled_value = resample(double(video_pc_traces(1:num_pcs,chunk_start:chunk_end))',...
                round(fps*mult_factor),round(resample_rate*mult_factor));
            video_pc_traces_resampled{lk}( frame_current+1:(frame_current+size(resampled_value,1)),:) = resampled_value;
            frames_total = frames_total+size(resampled_value,1);
            frame_current = frame_current+size(resampled_value,1);
            %              video_pc_traces_resampled{lk} = cat(1, video_pc_traces_resampled{lk},...
            %                  resample(double(video_pc_traces(1:num_pcs,chunk_start:chunk_end))',round(fps*mult_factor),round(resample_rate*mult_factor)));
        end
        video_pc_traces_resampled{lk} = video_pc_traces_resampled{lk}(1:frame_current,:);
        fclose(f)
        