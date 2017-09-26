run vl_compilenn ;
run vl_setupnn ;

% Download a pre-trained CNN from the web (needed once).
urlwrite(...
  'http://www.vlfeat.org/matconvnet/models/imagenet-vgg-f.mat', ...
  'imagenet-vgg-f.mat') ;

% Setup MatConvNet.
net = dagnn.DagNN.loadobj(load('imagenet-resnet-50-dag.mat')) ;
net.mode = 'test' ;

% Load a model and upgrade it to MatConvNet current version.
net = load('imagenet-vgg-f.mat') ;
net = vl_simplenn_tidy(net) ;

% Obtain and preprocess an image series
imageseries_folder = 'E:\Bence\Data\Motionanalysis_captures\20170228\Vicon6_longrun2\636238870803681649';
cameradirectory = 'E:\Bence\Data\Motionanalysis_captures\20170228\Vicon6_longrun2\CameraU';
read_mkv_file(imageseries_folder,cameradirectory );

%imds = imageDatastore(moviefiles(3).name,'ReadFcn',@mpfour_reader);

moviefiles = dir(strcat(cameradirectory,filesep));
videofeatures = [];


videofile = [];
for ll = 3:numel(moviefiles)
    imagefile = moviefiles(ll).name;
    temp = mpfour_reader(imagefile);
    videofile = cat(2,videofile,temp);
end

videomatrix = single(cat(4,videofile(:).cdata));
videomatrixnormalized = bsxfun(@minus,videomatrix,mean(videomatrix,4));

for kk = 1:size(videofile,2)
    
%im = imread('peppers.png') ;
im = squeeze(videomatrixnormalized(:,:,:,kk));
im_ = single(im) ; % note: 255 range
im_ = imresize(im_, net.meta.normalization.imageSize(1:2)) ;
im_ = im_ - net.meta.normalization.averageImage ;

% Run the CNN.
res = vl_simplenn(net, im_) ;

% % Show the classification result.
% scores = squeeze(gather(res(end).x)) ;
% [bestScore, best] = max(scores) ;
% figure(1) ; clf ; imagesc(im) ;
% title(sprintf('%s (%d), score %.3f',...
%    net.meta.classes.description{best}, best, bestScore)) ;

videofeatures =   cat(2,videofeatures,squeeze(gather(res(end-1).x))) ;

end
