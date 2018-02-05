
function neighborsCell = identifyNeighborsAuto(icmat, handles)

%%%%%%%%%%%%
% this code automatically sorts through to find all cell neighbors within a
% certain distance of the target (boundary to boundary). the output is a
% cell with the vector of the neighbor indices in the target cell's entry.
%
% Laurie Burns, Sept 2010.
%%%%%%%%%%%

global SpikeImageData

plottingOn = 0;

overlapradius = handles.OverlapRadius;

cellvec(1) = 1;
cellvec(2) = size(icmat,3);


%% look for overlap of IC and dilated IC
neighborsCell=cell(cellvec(2),1);
se = strel('disk',overlapradius,0);
for c = 1:cellvec(2)
    thisCellDilateCopy = repmat(imdilate(icmat(:,:,c),se),[1 1 size(icmat,3)]);
    res = icmat.*thisCellDilateCopy;
    res = squeeze(sum(sum(res,2),1));
    res = find(res>0);
    neighborsCell{c,1} = setdiff(res,c);
    
    if mod(c,50)==1
        fprintf('up to cell number %d \n',c)
    end
end

%% if want to plot it
if plottingOn
    figure;
    colormap(gray);axis image
    hold on
    for c=1:size(icmat,3)
        contour(gaussblur1(icmat(:,:,c),2),1)
        [x,y] = ait_centroid(icmat(:,:,c));
        text(x-2,y,num2str(c),'fontsize',10)
    end
    clear c x y
    
    for cellnum = cellvec(1):cellvec(2)
        % bold red the main cell
        [x,y] = ait_centroid(icmat(:,:,cellnum));
        hmain = text(x-2,y,num2str(cellnum),'fontsize',10,'fontweight','bold','color','r');
        %     handle_array = cell(10,1);
        handle_array = [];
        
        for d = 1:length(neighborsCell{cellnum,1})
            %         counter = counter + 1;
            c = neighborsCell{cellnum,1}(d);
            [x,y] = ait_centroid(icmat(:,:,c));
            h = text(x-2,y,num2str(c),'fontsize',10,...
                'fontweight','bold','color','b');
            handle_array{d,1} = h;
            %             neighborVector = cat(1,neighborVector,c);
            %         end
        end
        %     neighborsCell{cellnum,1} = neighborVector;
        pause()
        delete(hmain)
        for hnum = 1:length(handle_array)
            delete(handle_array{hnum,1})
        end
    end
end
