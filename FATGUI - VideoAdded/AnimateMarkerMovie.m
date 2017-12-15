function M = AnimateMarkerMovie(mocapstruct, frame_inds, hObject, eventdata, handles)

global CurrentFrame;

videofolder = mocapstruct.cameradirectory;
matchedindex = mocapstruct.matched_frames_aligned(frame_inds);

h1=handles.a1;
h2=handles.a2;

frame_last = 0;

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
% subplot(2,1,1) % ##
set(h1,'Color','k');

xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1));
yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2));
zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3));
line(h1,xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);

if CurrentFrame ==1
    
    axis(h1,'manual');
    set(h1,'Color','k');
    grid(h1,'on');
    
    set(h1,'Xcolor',[1 1 1 ]);
    set(h1,'Ycolor',[1 1 1]);
    set(h1,'Zcolor',[1 1 1]);

    zlim(h1,[-100 150]);
    xlim(h1,[-300 300]);
    ylim(h1,[-300 300]);
    set(h1,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]);
    view(h1,[-22, 12]);
  

    
    imagesc(h2,[]);
    set(h2,'Color',[0.94 0.94 0.94]);
    set(h2,'Xcolor',[0.94 0.94 0.94]);
    set(h2,'Ycolor',[0.94 0.94 0.94]);
    set(h2,'Zcolor',[0.94 0.94 0.94]);
  
    set(h2,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]);
    set(h2,'XTick',[],'YTick',[],'ZTick',[]);
    

end

oldvidobj = [];
oldfile = 'oldfile.mp4';
oldframe = 0;
fighandle = [];
base_time = datenum(mocapstruct.mocapfiletimes{1});

for lk = reshape(frame_inds,1,[])
    
    cla(h1);
    %cla(h2);
    
    %% Plot the marker movie
    set(handles.t4,'String',num2str(lk));
    ind_to_plot = lk;
    %subplot(2,1,1) % ##
    
    %% Plot markers that are tracked in the frame
    set(h1,'Nextplot','ReplaceChildren');
    set(h2,'Nextplot','ReplaceChildren');
    handles_here = cell(1,numel(mocapstruct.markernames));
    for jj = 1:numel(mocapstruct.markernames)
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
            if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
              
                xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                handles_here{jj} = line(h1,xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',8);

                hold(h1,'on');
                marker_plot(jj) = 1;
            else
                marker_plot(jj) = 0;
            end
        end
    end
    
    if marker_plot(5) ==0 || marker_plot(4) == 0
        marker_plot(:) = 0;
    end
    
    %% plot the links between markers
    for mm = 1:numel(mocapstruct.links)
        if (ismember(mocapstruct.links{mm}(1),1:numel(mocapstruct.markernames)) && ismember(mocapstruct.links{mm}(2),1:numel(mocapstruct.markernames)))
            if (marker_plot(mocapstruct.links{mm}(1)) == 1 && marker_plot(mocapstruct.links{mm}(2)) == 1)

                xx = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
                yy = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
                zz = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
                
                line(h1,xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',3);  
            end
        end
    end
    
    title(h1,strcat('  Frame: ' ,datestr(base_time+lk./(mocapstruct.fps*60*60*24))),'Color','w');
    
    
    drawnow;
    hold(h1,'off');
    
    frame_last = lk;
    
    %% plot the video movie
    videoframe = matchedindex(find(frame_inds==lk));
    framenumber = mod(videoframe,3500);
    
    if (videoframe>0 && oldframe ~=videoframe)
        
        videofile = strcat(videofolder,num2str(videoframe-framenumber),'.mp4');
        
        % subplot(2,1,2); % ##
        [imagesout,oldvidobj,oldfile,fighandle] = mpfour_reader_singleframe(videofile,framenumber,oldvidobj,oldfile,oldframe,fighandle);
        imagesc(h2,imagesout(1).cdata);
        oldframe = framenumber;
        set(h2,'xtick',[],'ytick',[]);
        
    end
    
    M(find(frame_inds == lk)) =  getframe(h2);
    
end

set(h1,'Nextplot','add');
set(h2,'Nextplot','add');

end