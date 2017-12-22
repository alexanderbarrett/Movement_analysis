function make_cnn_features(videodirectory_in)
% Load a model and upgrade it to MatConvNet current version.
% net = load('imagenet-vgg-f.mat') ;
% net = vl_simplenn_tidy(net) ;

% Long folder
%imageseries_folder = 'E:\Bence\Data\Motionanalysis_captures\20170228\Vicon6_longrun2\636238870803681649';
%cameradirectory_U = 'E:\Bence\Data\Motionanalysis_captures\20170228\Vicon6_longrun2\CameraU';

%1 hr long recording
% imageseries_folder_U = 'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\636235234683982964';
% cameradirectory_U = 'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU';
%
% imageseries_folder_L = 'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\636235234727205764';
% cameradirectory_L = 'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL';
camera_base_agg = [];
%camera_base_agg{1} = 'E:\Bence\Data\Motionanalysis_captures\20170224\Toydog_long1';
%camera_base_agg{2} = 'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1';
%camera_base_agg{9} = 'E:\Bence\Data\Motionanalysis_captures\20170228\Vicon6_longrun2';
%camera_base_agg{3} = 'E:\Bence\Data\Motionanalysis_captures\20170301\Vicon6_postanesthesia_ketamine1';
%camera_base_agg{5} = 'E:\Bence\Data\Motionanalysis_captures\20170302\vicon6_twotap_mid3';
%  camera_base_agg{6} = 'E:\Bence\Data\Motionanalysis_captures\20170227\Vicon6_pretask1';
%  camera_base_agg{7} = 'E:\Bence\Data\Motionanalysis_captures\20170227\Vicon6_task1';
%  camera_base_agg{8} = 'E:\Bence\Data\Motionanalysis_captures\20170227\Vicon6_novelobj1';
 camera_base_agg{10} = videodirectory_in;
 %camera_base_agg{11} = 'E:\Bence\Data\Motionanalysis_captures\20170309\toydog_no_markers';


createmasks(camera_base_agg,3,0)
overwrite_features = 1;
do_mask = 1;
run vl_compilenn ;
run vl_setupnn ;
vl_testnn('gpu', true)

% Download a pre-trained CNN from the web (needed once).
% urlwrite(...
%   'http://www.vlfeat.org/matconvnet/models/imagenet-vgg-f.mat', ...
%   'imagenet-vgg-f.mat') ;

% Setup MatConvNet.
net = dagnn.DagNN.loadobj(load('imagenet-resnet-50-dag.mat')) ;
net.mode = 'test' ;
net.removeLayer({'prob'});

do_conv = 1;
do_pca_layer = 0;
layer_num = 90; %middle
load_coefficients = 0;
%layer_num = 6; %start
%output layer means subtracting just 1, first layer go donw to 5
if (do_conv)
    for k=numel(net.layers):-1:layer_num
        net.removeLayer(net.layers(k).name);
    end
end

do_gpu = 0;
if (do_gpu)
    net.move('gpu');
end

base_agg_list = [10];

num_features_to_save = 1000;


for klk = base_agg_list
    %% turn this into a separate modular file
    for cameratype_ind = 1:3
        suffix = {'U','L','R'};
        cameradirectory = strcat(camera_base_agg{klk},filesep,'Camera',suffix{cameratype_ind});
        %if extension exist, lists subdirectories
        if exist(cameradirectory,'dir')
            listing = dir(cameradirectory);
            cameradirectory = strcat(cameradirectory,filesep,listing(3).name,filesep);
            
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
        
        maskfilename = strcat(cameradirectory,filesep,'mask_vals.mat');
        load(maskfilename);
        
        filelist = 1:numel(moviefiles);
        for ll = 1%filelist%
            movie_name_here = strcat(cameradirectory,filesep,num2str(vals(ll)),'.mkv');
            fprintf('For videofile %s loading in images \n',movie_name_here)
            tic
            videofile = [];
            %  imagefile = strcat(cameradirectory,filesep,movie_name_here);
            temp = mpfour_reader(movie_name_here,[]);
            videofile = cat(2,videofile,temp);
            
            videomatrix = single(cat(4,videofile(:).cdata));
            if (ll ==filelist(1))
                mean_videomatrix = mean(videomatrix,4);
            end
            videomatrixnormalized = bsxfun(@minus,videomatrix,mean_videomatrix);
            mask_position = round(mask_position);
            videomatrixnormalized = videomatrixnormalized(mask_position(2):(mask_position(2)+mask_position(4)),...
                mask_position(1):(mask_position(1)+mask_position(3)),:,:);
            toc
            
            fprintf('For videofile %f starting feature creation \n',ll)
            tic
            scores_intermediate = [];
            for kk = 1:size(videofile,2)
                
                %im = imread('peppers.png') ;
                im = squeeze(videomatrixnormalized(:,:,:,kk));
                im_ = single(im) ; % note: 255 range
                im_ = imresize(im_, net.meta.normalization.imageSize(1:2)) ;
                im_ = im_ - net.meta.normalization.averageImage ;
                % Run the CNN.
                %res = vl_simplenn(net, im_) ;
                if (do_gpu)
                    im_gpu = gpuArray(im_);
                    net.eval({'data', im_gpu}) ;
                else
                    net.eval({'data', im}) ;
                end
                % % Show the classification result.
                % scores = squeeze(gather(res(end).x)) ;
                % [bestScore, best] = max(scores) ;
                % figure(1) ; clf ; imagesc(im) ;
                % title(sprintf('%s (%d), score %.3f',...
                %    net.meta.classes.description{best}, best, bestScore)) ;
                %  scores = net.vars(end).value ;
                % test =squeeze(gather(scores));
                %scores = net.vars(end-1).value ;
                
                if (kk ==1)
                    score_val = gather(net.vars(end).value);
                    scores_intermediate = zeros(size(score_val,1),size(score_val,2),size(score_val,3),size(videofile,2));
                end
                scores_intermediate(:,:,:,kk) =   squeeze(gather(net.vars(end).value)) ;
                %                 if (do_conv)
                %
                %                 else
                %                     scores_intermediate =   cat(2,scores_intermediate,squeeze(gather(net.vars(end).value))) ;
                %                 end
                %videofeatures =   cat(2,videofeatures,squeeze(gather(res(end-1).x))) ;
            end
            toc
            fprintf('running pca and sving \n')
            
            
            if (do_pca_layer)
                tic
                test = reshape(scores_intermediate,size(scores_intermediate,1)*size(scores_intermediate,2)*size(scores_intermediate,3),size(scores_intermediate,4));
                
                                features_base = strcat(cameradirectory,filesep,'videofeatures_convlayer_pca_',num2str(layer_num),'_');

                if (ll == filelist(1))
                    [pca_coeff,~,~] = pca(test');
                    save(strcat(cameradirectory,filesep,features_base,'pca_coefficients'),'pca_coeff');
                end
                reconstructed_score = test'*pca_coeff;
                feature_matrix = reconstructed_score(:,1:num_features_to_save);
                featuresfile = strcat(features_base,num2str(ll),'_',...
                    suffix{cameratype_ind},'.mat');
                if (~exist(featuresfile,'file') || overwrite_features)
                  %  save(featuresfile,'feature_matrix','-v7.3');
                end
                toc
                
            else
                % if taking only the last layer don't worry too much about
                % the size
                   
                featuresfile = strcat(cameradirectory,filesep,'videofeatures_outputlayer',num2str(ll),'_',...
                    suffix{cameratype_ind},'.mat');
                if (~exist(featuresfile,'file') || overwrite_features)
                   save(featuresfile,'scores_intermediate','-v7.3');
                end
                toc
                
              %  videofeatures =   cat(4,videofeatures,scores_intermediate) ;
                
                
            end
            % save to separate files

        end
        
%         if (~do_pca_layer)
%             if (~exist(featuresfile,'file') || overwrite_features)
%                 save(featuresfile,'videofeatures');
%             end
%             
%         end
        
        
        %         figure(3000+klk+10*cameratype_ind)
        %         imagesc(squeeze(videofeatures)')
        %         featuresfile = strcat(cameradirectory,filesep,'videofeatures',suffix{cameratype_ind},'.mat');
        end
    end
end
end