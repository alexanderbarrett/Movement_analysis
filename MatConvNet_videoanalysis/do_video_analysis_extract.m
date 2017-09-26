
if (do_video_analysis)
    %loop over different possibilities fro number ofinitial frames missing
    %to see which improves fit
    video_fps = 60;
    %video_missingframes = floor((50/60)*fps); %note this is in terms of the normal traces
    video_pc_traces_resampled = cell(1,numel(videofeatures_directory));
    inds_to_fit = cell(1,numel(videofeatures_directory));
    video_lengths = zeros(1,numel(videofeatures_directory));
    
    
    
    %% load in the video features or create them if not present
    for lk = 1:numel(videofeatures_directory)
        %% get the .times file
        %% load in the video frames
        times_files = dir(strcat(videofeatures_directory{lk},filesep,'*.times'));
        
        
        
        
        %% if using the last layer or intermediate layesr
        %             if (~do_conv_features)
        %                            videofile_here = videofeatures_files{lk};
        %          videofeatures = (load(videofile_here));
        %             videofeatures_here = squeeze(gather(videofeatures.videofeatures));
        %
        % else
        %  tag = 'convlayer*';
        % videostart = 'videofeatures';
        %% need to standardize name
        videofilenames = dir(strcat(videofeatures_directory{lk},filesep,feature_tag));
        % videofilenames = dir(strcat(videofeatures_directory{lk},filesep,videostart,'_',tag));
        fprintf('Loading video files for camera %f \n',lk)
        videofeatures_agg = [];
        load_inds = zeros(1,numel(videofilenames));
        for ll= 1:numel(videofilenames)
            % thanks Jan Simon!
            Str = videofilenames(ll).name;
            Key   = strrep(feature_tag,'*','');
            Index = strfind(Str, Key);
            load_inds(ll) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
        end
        [~,ind_to_load] = sort(load_inds,'ascend');
        
        for ll = ind_to_load
            fprintf('loading file %s \n',videofilenames(ll).name);
            test = load(strcat(videofeatures_directory{lk},filesep,videofilenames(ll).name));
            
            fieldname_here = fieldnames(test);
            [~,fieldagg] = max(size(getfield(test,fieldname_here{1})));
            videofeatures_agg = cat(fieldagg,videofeatures_agg, getfield(test,fieldname_here{1}));
        end
        %% take precautions here and above to definte the correct dimension
        [~,fieldmax] = max(size(squeeze(videofeatures_agg)));
        if (fieldmax == 1)
            videofeatures_agg = videofeatures_agg';
        end
        
        videofeatures_here = squeeze(videofeatures_agg);
        fprintf('finished loading video files \n')
        
        original_length = size(videofeatures_here,2);
        %for screw up
        frames_to_analyze = min(original_length,70000);
        
        videofeatures_here = videofeatures_here(:,1: frames_to_analyze);
        
        %get the PCs of the features
        [coeff,score,latent,tsquared,explained] = pca(squeeze(videofeatures_here)');
        video_pc_traces = score';
        %video_pc_traces = videofeatures_here;
        %resample the pcs
        num_pcs = 150;
        
        
        %% resample at the variable framerate
        %true_video_framerate = fps*(original_length+450+video_missingframes*(video_fps./fps))./size(markers_preproc.(marker_names{1}),1);
        true_video_framerate = 59.82;
        % video_pc_traces_resampled{lk} = resample(double(video_pc_traces(1:num_pcs,:))',round(fps*100),round(true_video_framerate*100));
        
        f = fopen(strcat(videofeatures_directory{lk},filesep,times_files.name));
        float1 = fread(f,[1,100000000],'uint64');
        frame_number = numel(float1);
        % 1000.*(frame_number./(float1(end)-float1(1)));
        frame_chunk = 5000;
        
        %preallocate
        video_pc_traces_resampled{lk} = zeros( floor(frames_to_analyze*245./true_video_framerate*1.1),num_pcs);
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
        
        %% subtract the missing frames from the mocap trace
        video_length = size(video_pc_traces_resampled{lk},1);
        video_lengths(lk) = video_length;
        
        %use the PCs to learn a mapping to marker positions
        figure(2300+lk)
        imagesc(squeeze(video_pc_traces_resampled{lk})')
        
    end
    
    
    %% loop over different 'microoffsets'
    missingframes_here = [1];
    %agg_mse = cell(1,numel(missingframes_here));
    do_plots_regression = 0;
    for video_missingframes = missingframes_here;
        fprintf('starting missing frames %f \n',video_missingframes);
        for lk = 1:numel(videofeatures_directory)
            
            inds_to_fit{lk} = (video_missingframes:(video_lengths(lk) +video_missingframes-1));
        end
        
        %% obtain regressors
        if (numel(video_lengths)>1)
            [min_length_val,min_ind] = min(video_lengths);
            ind_to_fit = inds_to_fit{min_ind};
            regression_input = video_pc_traces_resampled{1}(1:min_length_val,:);
            for lk = 2:numel(video_lengths)
                regression_input = cat(2, regression_input,...
                    video_pc_traces_resampled{lk}(1:min_length_val,:));
            end
            
        else
            ind_to_fit = inds_to_fit{1};
            regression_input = video_pc_traces_resampled{1};
        end
        
        %% save regression input
        
        % choose features to cluster over
        
        %do individual regressions for each relative marker position
        %Mdl = fitrsvm(video_pc_traces_resampled, delta_markers_reshaped(ind_to_fit,1),'Kernelfunction','rbf' );
        marker_scan = 1:numel(marker_names);
        num_markers = numel(marker_scan);
        num_poses = 1;
        do_regression_markers = 1;
        
        %% create better output features
        agg_features = [];
        % feature_mapping = cell(1,numel(marker_scan)*3);
        for ll = marker_scan
            agg_features = cat(1,agg_features,markers_preproc.(marker_names{ll})(:,:)');
        end
        [feature_coeff,feature_score,latent,tsquared,explained] = pca(agg_features','Centered','off');
        num_features_train = size(feature_score,2);
        %% break up training into different poses, based on ind to fit
        pose_indicies = cell(1,num_poses);
        headpos = markers_preproc.(marker_names{1})(ind_to_fit,:);
        headpos = bsxfun(@minus,headpos,median(headpos,1));
        headangle = atan2d(headpos(:,1),headpos(:,2));
        
        figure(1233)
        subplot(2,1,1)
        plot(headpos(:,1),headpos(:,2))
        subplot(2,1,2)
        plot(headangle)
        
        angle_bins = -180:360./(num_poses):180;
        
        headangle_binarized = zeros(size(headangle));
        
        for lk = 1:num_poses
            
            pose_indicies{lk} = sort(intersect(find(headangle>angle_bins(lk)),...
                find(headangle<=angle_bins(lk+1))),'ascend');
            headangle_binarized(intersect(find(headangle>angle_bins(lk)),...
                find(headangle<=angle_bins(lk+1)))) = lk;
        end
        num_outputs = size(markers_preproc.(marker_names{1}),2);
        
        if (do_regression_markers)
            Mdl = cell(num_markers,num_outputs,num_poses);
            fitinfo = cell(num_markers,num_outputs,num_poses);
            model_mse = cell(2,num_markers,num_outputs,num_poses);
        else
            Mdl = cell(num_features_train,1,num_poses);
            fitinfo = cell(num_features_train,1,num_poses);
            model_mse = cell(2,num_features_train,1,num_poses);
        end
        %don't overwrite these
        regression_input_original = regression_input;
        ind_to_fit_original = ind_to_fit;
        
        movies = cell(2,2,num_poses); %for real/regressed one for input, one for output
        
        for lk = 1:num_poses
            fprintf('for pose %f \n',lk);
            regression_input = regression_input_original;
            ind_to_fit = ind_to_fit_original;
            % restrict to a given pose
            regression_input = regression_input(pose_indicies{lk},:);
            ind_to_fit = ind_to_fit(pose_indicies{lk});
            
            input_length = size(regression_input,1);
            training_data = [1:floor(input_length *0.4),(floor(input_length*0.5)):input_length];
            held_out_data = (1+floor(input_length *0.4)):(floor(input_length*0.5));
            kfold_training_data= 1:floor(input_length);
            
            regression_output = markers_preproc.(marker_names{1})(ind_to_fit,:);
            kfold_predicted_score = zeros(numel(kfold_training_data) ,num_features_train);
            traintest_predicted_score = cell(1,2);
            
            do_kfold = 1;
            
            if (~do_regression_markers)
                %score*feature_coeff'
                for jj = 1:num_features_train
                    fprintf('training for feature %f \n',jj)
                    regression_output = feature_score(ind_to_fit,jj);
                    
                    if (do_kfold)
                        Mdl{jj,1,lk} = fitrlinear(regression_input(kfold_training_data,:), squeeze(regression_output(kfold_training_data)),...
                            'Learner','leastsquares','Solver','sparsa','GradientTolerance',1e-6,'BetaTolerance',1e-6,'Crossval','on' ...
                            );
                        model_mse(1,jj,1,lk) = kfoldLoss(Mdl{jj,1,lk});
                        % if (do_plots_regression)
                        figure(2000+jj+10000*lk)
                        %   subplot(2,num_outputs,1)
                        reduce_plot(Mdl{jj,1,lk}.kfoldPredict ,'+r')
                        hold on
                        reduce_plot(regression_output(kfold_training_data),'+k')
                        hold off
                        kfold_predicted_score(:,jj,1) = Mdl{jj,1,lk}.kfoldPredict;
                        
                    else
                        %                     Mdl{jj,1,lk} = fitrlinear(regression_input(training_data ,:), squeeze(regression_output(training_data )),...
                        %                         'Learner','leastsquares','Solver','sparsa','GradientTolerance',1e-6,'BetaTolerance',1e-6...
                        %                         );
                        %
                        Mdl{jj,1,lk} = fitrtree( regression_input(training_data ,:), squeeze(regression_output(training_data )));
                        
                        model_mse(1,jj,1,lk) = sqrt(Mdl{jj,1,lk}.loss(regression_input(training_data,:),regression_output(training_data,1),'LossFun'   ,'mse'));
                        model_mse(2,jj,1,lk) = sqrt(Mdl{jj,1,lk}.loss(regression_input(held_out_data,:),regression_output(held_out_data,1),'LossFun'   ,'mse'));
                        prediction{ll} = Mdl{jj,1,lk}.predict(regression_input(training_data,:));
                        
                        traintest_predicted_score{1} = cat(2, traintest_predicted_score{1},Mdl{jj,1,lk}.predict(regression_input(training_data,:)));
                        traintest_predicted_score{2} = cat(2, traintest_predicted_score{2},Mdl{jj,1,lk}.predict(regression_input(held_out_data,:)));
                        
                        %kFoldPredict(regression_input(training_data,:));
                        if (do_plots_regression)
                            figure(1000+jj+10000*lk)
                            subplot(2,1,1)
                            % reduce_plot(test ,'+r')
                            
                            reduce_plot(prediction{ll} ,'+r')
                            hold on
                            reduce_plot(regression_output(training_data),'+k')
                            hold off
                            
                            prediction{ll} = Mdl{jj,1,lk}.predict(regression_input(held_out_data,:));
                            
                            subplot(2,1,2)
                            reduce_plot(prediction{ll},'+r')
                            hold on
                            reduce_plot(regression_output(held_out_data),'+k')
                            hold off
                        end
                    end
                end
                
                
                reconstructed_pose_true  = feature_score*feature_coeff'; %agg_features';%
                
                for mm =1:2
                    
                    if (do_kfold)
                        reconstructed_pose = kfold_predicted_score*feature_coeff';
                    else
                        reconstructed_pose = traintest_predicted_score{mm}*feature_coeff';
                    end
                    
                    
                    if (mm ==1)
                        data_use_here = training_data;
                    else
                        data_use_here = held_out_data;
                    end
                    
                    %                     for jj = marker_scan
                    %                         markers_preproc_true.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(data_use_here),:);
                    %                         test = zeros(numel(data_use_here),num_outputs);
                    %                         for ll = 1:num_outputs
                    %                             test(:,ll) = Mdl{jj,ll,lk}.predict(regression_input(data_use_here,:));
                    %                         end
                    %                         markers_preproc_predicted.(marker_names{jj}) = test;
                    %                         %markers_preproc.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(held_out_data));
                    %                     end
                    %
                    reconstructed_markers = markers_preproc;
                    reconstructed_markers_true = markers_preproc;
                    
                    for zz =marker_scan
                        for running_ind = 1:3
                            reconstructed_markers.(marker_names{zz})(1:numel(data_use_here),(running_ind)) =reconstructed_pose(1:numel(data_use_here),3*(zz-1)+running_ind);
                            reconstructed_markers_true.(marker_names{zz})(1:numel(data_use_here),(running_ind)) =reconstructed_pose_true(ind_to_fit(data_use_here),3*(zz-1)+running_ind);
                        end
                    end
                    
                    
                    matlab_fr = 1;
                    frames_to_use = 25000;
                    trace_to_use = 1:numel(data_use_here);
                    frame_inds = trace_to_use(1:matlab_fr:min(frames_to_use,numel(trace_to_use)))';
                    save_movie = 1;
                    figure(370)
                    %initialize movies
                    movies{1,mm,lk}(1) = getframe(gcf);
                    movies{2,mm,lk}(1) = getframe(gcf);
                    
                    movies{1,mm,lk} = animate_markers(reconstructed_markers,frame_inds,marker_names(marker_scan),markercolor,links,movies{1,mm,lk},save_movie);
                    movies{2,mm,lk} = animate_markers(reconstructed_markers_true,frame_inds,marker_names(marker_scan),markercolor,links,movies{2,mm,lk},save_movie);
                    
                end
                
                
            else
                % specify_num_workers(8);
                %parpool('local',4)
                for jj =9%:numel(marker_names)
                    fprintf('for marker %f \n',jj);
                    %regression_input = cat(2,video_pc_traces_resampled(:));
                    regression_output = markers_preproc.(marker_names{jj})(ind_to_fit,:);
                    %regression_input = cat(2,regression_input, reshape(marker_position(setxor(1:size(marker_position,1),jj),ind_to_fit,:),numel(ind_to_fit),[]) );
                    prediction = cell(1,num_outputs);
                    
                    
                    for ll = 1%2:num_outputs
                        hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus');
                        Lambda = logspace(-5,-1,15);
                        % [Mdl{jj,ll},fitinfo{jj,ll}] = fitrlinear(regression_input(training_data,:), squeeze(regression_output(training_data,ll)),...
                        %    'Lambda',Lambda, 'Learner','leastsquares','Solver','sparsa','Regularization','lasso','GradientTolerance',1e-6,'BetaTolerance',1e-6 ...
                        %      );
                        %                         [Mdl{jj,ll,lk},fitinfo{jj,ll,lk}] = fitrlinear(regression_input(training_data,:), squeeze(regression_output(training_data,ll)),...
                        %                             'Learner','leastsquares','Solver','sparsa','GradientTolerance',1e-6,'BetaTolerance',1e-6 ...
                        %                             );
                        fprintf('starting cross validated regression for output: %f marker %f \n',ll,jj);
                        if (do_kfold)
                            Mdl{jj,ll,lk} = fitrtree( [regression_input(kfold_training_data(1:100000) ,:)], squeeze(regression_output(kfold_training_data(1:100000),ll )),...
                                'KFold',5 );
                            
                            Mdl_nocv = fitrtree( [regression_input(kfold_training_data(1:100000) ,:)], squeeze(regression_output(kfold_training_data(1:100000),ll )),...
                                'MaxNumSplits',20,'MinLeafSize',10,'KFold',5);
                            
                            loss( Mdl_nocv,regression_input(kfold_training_data(1:100000) ,:),squeeze(regression_output(kfold_training_data(1:100000),ll )),'LossFun','mse')
                            loss(Mdl_nocv,regression_input(kfold_training_data(100000:110000) ,:),squeeze(regression_output(kfold_training_data(100000:110000),ll )),'LossFun','mse')
                            
                            Md
                            %,'MaxNumSplits',7 MinLeafSize, or MinParentSize
                            %Mdl =
                            %
                            Mdl_test = TreeBagger(100,regression_input(kfold_training_data([1:100000 120000:250000]) ,:),...
                                squeeze(regression_output(kfold_training_data([1:100000 120000:250000]),ll )),'Method','regression','OOBPrediction','on');
                            
                            Mdl_test.oobError()
                            Mdl_test.error(regression_input(kfold_training_data([1:100000 120000:250000]) ,:),...
                                squeeze(regression_output(kfold_training_data([1:100000 120000:250000]),ll )));
                            
                            Mdl_test.error(regression_input(kfold_training_data(100000:110000) ,:),...
                                squeeze(regression_output(kfold_training_data(100000:110000),ll )));
                            
                            testval = Mdl_test.predict(regression_input(kfold_training_data(1:100000) ,:));
                            sqrt(mean( ( squeeze(regression_output(kfold_training_data(1:100000),ll ))-testval).^2))
                            
                            testval = Mdl_test.predict(regression_input(kfold_training_data(105000:115000) ,:));
                            sqrt(mean( ( squeeze(regression_output(kfold_training_data(105000:115000),ll ))-testval).^2))
                            %                                                   Mdl_test = fitrtree( [headangle_binarized regression_input(kfold_training_data ,:)], squeeze(regression_output(kfold_training_data )),...
                            %                                                       'KFold',5,'CategoricalPredictors',1 );
                            
                            model_mse{1,jj,ll,lk} = sqrt(Mdl{jj,ll,lk}.kfoldLoss('lossfun'   ,'mse'));
                            %model_mse(2,jj,ll,lk) = sqrt(Mdl{jj,ll,lk}.kfoldLoss('lossFun'   ,'mse'));
                            prediction{ll} = Mdl{jj,ll,lk}.kfoldPredict();
                            
                        else
                            Mdl{jj,ll,lk} = fitrtree( regression_input(training_data ,:), squeeze(regression_output(training_data )));
                            %model_mse{1,jj,ll,lk} = sqrt(Mdl{jj,ll,lk}.loss(regression_input(training_data,:),regression_output(training_data,ll),'LossFun'   ,'mse'));
                            %  model_mse{2,jj,ll,lk} = sqrt(Mdl{jj,ll,lk}.loss(regression_input(held_out_data,:),regression_output(held_out_data,ll),'LossFun'   ,'mse'));
                            prediction{ll} = Mdl{jj,ll,lk}.predict(regression_input(training_data,:));
                        end
                        %
                        %
                        % [Mdl{jj,ll}] = fitrsvm(regression_input(1:10:30000,:), squeeze(regression_output(1:10:30000,ll)),...
                        %     'KernelFunction','rbf','CrossVal' );
                        
                        
                        
                        %kFoldPredict(regression_input(training_data,:));
                        if (do_plots_regression)
                            if (~do_kfold)
                                figure(1000+jj+10000*lk+20*video_missingframes)
                                subplot(2,num_outputs,ll)
                                % reduce_plot(test ,'+r')
                                
                                reduce_plot(prediction{ll} ,'+r');
                                hold on
                                reduce_plot(regression_output(training_data,ll),'+k');
                                hold off
                                
                                
                                prediction{ll} = Mdl{jj,ll,lk}.predict(regression_input(held_out_data,:));
                                
                                subplot(2,num_outputs,num_outputs+ll)
                                reduce_plot(prediction{ll},'+r');
                                hold on
                                reduce_plot(regression_output(held_out_data,ll),'+k');
                                hold off
                            else
                                figure(1000+jj+10000*lk+20*video_missingframes)
                                subplot(1,num_outputs,ll)
                                % reduce_plot(test ,'+r')
                                
                                reduce_plot(prediction{ll} ,'+r');
                                hold on
                                reduce_plot(medfilt1(prediction{ll},5) ,'+b');
                                reduce_plot(regression_output(kfold_training_data,ll),'+k');
                                hold off
                                
                                
                                
                            end
                        end
                        
                    end
                end
            end
            %video_missingframes = missingframes_here;
            %agg_mse{find(video_missingframes == missingframes_here)} = model_mse;
            
            %% animate prediction model
            do_animation = 1;
            if do_animation
                
                mm_scan = 1:2;
                if (do_kfold)
                    mm_scan = 1;
                end
                for mm = mm_scan
                    if (~do_kfold)
                        if (mm ==1)
                            data_use_here = training_data;
                        else
                            data_use_here = held_out_data;
                        end
                    else
                        data_use_here = kfold_training_data;
                    end
                    
                    
                    markers_preproc_true = markers_preproc;
                    markers_preproc_predicted = markers_preproc;
                    
                    for jj = marker_scan
                        markers_preproc_true.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(data_use_here),:);
                        test = zeros(numel(data_use_here),num_outputs);
                        for ll = 1:num_outputs
                            if (~do_kfold)
                                test(:,ll) = Mdl{jj,ll,lk}.predict(regression_input(data_use_here,:));
                            else
                                test(:,ll) = Mdl{jj,ll,lk}.kfoldPredict();
                            end
                        end
                        markers_preproc_predicted.(marker_names{jj}) = test;
                        %markers_preproc.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(held_out_data));
                        
                        
                    end
                    %markers_preproc.(marker_names{1})
                    %movies{1,1,lk}
                    %movies{1,1,lk} = movie;
                    matlab_fr = 10;
                    trace_to_use = 1:size(markers_preproc_true.(marker_names{jj}),1);
                    frames_to_use = 30000;%numel(trace_to_use);
                    frame_inds = trace_to_use(1:matlab_fr:min(frames_to_use,numel(trace_to_use)))';
                    save_movie = 1;
                    figure(370)
                    %initialize movies
                    movies{1,mm,lk}(1) = getframe(gcf);
                    movies{2,mm,lk}(1) = getframe(gcf);
                    
                    movies{2,mm,lk} = animate_markers(markers_preproc_predicted,frame_inds,marker_names(marker_scan),markercolor,links,movies{2,mm,lk},save_movie);
                    movies{1,mm,lk} = animate_markers(markers_preproc_true,frame_inds,marker_names(marker_scan),markercolor,links,movies{1,mm,lk},save_movie);
                end
                
                
                
                save_movie_animation = 1;
                if (save_movie_animation)
                    save_tags{1,1} = 'true_training_2.mp4';
                    save_tags{2,1} = 'predicted_training_2.mp4';
                    save_tags{1,2} = 'true_testing.mp4';
                    save_tags{2,2} = 'predicted_testing.mp4';
                    
                    savedirectory_subcluster =strcat(savedirectory, save_tag,num2str(lk),filesep);
                    if (~exist(savedirectory_subcluster,'dir'))
                        mkdir(savedirectory_subcluster)
                    end
                    
                    for pp = 1:2
                        for mm = 1
                            if (numel(movies{pp,mm,lk}))
                                v = VideoWriter(strcat(savedirectory_subcluster,save_tags{pp,mm}),'MPEG-4');
                                open(v)
                                writeVideo(v,movies{pp,mm,lk})
                                close(v)
                            end
                        end
                    end
                end
            end
            fprintf('for pose %f model mse for missing frames  %f \n',lk,video_missingframes);
            %    model_mse(1,:,:,lk)
            %  model_mse(2,:,:,lk)
            
        end
        save_tag_model = strcat('Tree_regression_model_kfold',num2str(do_kfold),'_pcs_',num2str(num_pcs),'_poses_',num2str(num_poses),'_markers',num2str(numel(marker_scan)),'.mat');
        save_tag_mse = strcat('Tree_regression_model_kfold',num2str(do_kfold),'_pcs_',num2str(num_pcs),'_poses_',num2str(num_poses),'_markers',num2str(numel(marker_scan)),'mse.mat');
        save_directory_here = strcat(savedirectory,save_tag,filesep);
        if (~exist(save_directory_here,'dir'))
            mkdir(save_directory_here);
        end
        savedirectory_model =strcat(save_directory_here, save_tag_model);
        savedirectory_mse =strcat(save_directory_here, save_tag_mse);
        save(savedirectory_mse,'model_mse')
        save(savedirectory_model,'Mdl','-v7.3')
        for lk = 1:num_poses
            fprintf('for pose %f model mse for missing frames  %f \n',lk,video_missingframes);
            % model_mse(1,:,:,lk)
            %model_mse(2,:,:,lk)
        end
        
    end
end