function [mergedImage, cellColors] = mergeImages(cellImages,varargin)

options.colorType='random';
options=getOptions(options,varargin);


ColorMat=get(0,'DefaultAxesColorOrder');
imgSize=size(cellImages(:,:,1));
nImages=size(cellImages,3);
cellColors=zeros(nImages,3);

mergedImage=zeros([imgSize(1:2) 3]);
for i=1:nImages
    switch options.colorType
        case 'fixed'
            if i==1
                mergedImage=zeros([imgSize(1:2) 3]);
            end
            
            thisImage=double(cellImages(:,:,i));
            thisImage(thisImage<=0)=0;
            thisImage=(thisImage-min(min(thisImage)))/(max(max(thisImage))-min(min(thisImage)));
            
            IndexColor=mod(i,size(ColorMat,1))+1;
            colorLocal=ColorMat(IndexColor,:);
            cellColors(i,:)=colorLocal;
            
            colorLocalHSV=rgb2hsv(colorLocal);
            ToMax(:,:,1)=mergedImage(:,:,2);
            ToMax(:,:,2)=thisImage;
            [mergedImage(:,:,2),IndexSat]=max(ToMax,[],3);
            OldHue=mergedImage(:,:,1);
            OldHue(IndexSat==2)=colorLocalHSV(1);
            mergedImage(:,:,1)=OldHue;
            mergedImage(:,:,3)=mergedImage(:,:,2);
            
        case 'random'
            if i==1
                mergedImage=zeros([imgSize(1:2) 3]);
            end
            
            thisImage=double(cellImages(:,:,i));
            thisImage(thisImage<=0)=0;
            thisImage=(thisImage-min(min(thisImage)))/(max(max(thisImage))-min(min(thisImage)));
            
            colorLocal=rand([1 3]);
            cellColors(i,:)=colorLocal;
            
            colorLocalHSV=rgb2hsv(colorLocal);
            ToMax(:,:,1)=mergedImage(:,:,2);
            ToMax(:,:,2)=thisImage;
            [mergedImage(:,:,2),IndexSat]=max(ToMax,[],3);
            OldHue=mergedImage(:,:,1);
            OldHue(IndexSat==2)=colorLocalHSV(1);
            mergedImage(:,:,1)=OldHue;
            mergedImage(:,:,3)=mergedImage(:,:,2);
    end
end

mergedImage=hsv2rgb(mergedImage);
figure; imagesc(mergedImage)