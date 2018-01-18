function AnimateMarkerMovie(mocapstruct, frame_inds, hObject, eventdata, handles)

% Get frame offset value - base = 300,000
base = str2double(get(handles.eBase,'string'));

if isnan(base)
    base = 0;
end

frame_inds = frame_inds + base;

global videoPathValid;
if ~videoPathValid
    
    set(handles.tg1,'Value', 0);
    set(handles.tg2,'Value', 0);
    set(handles.tg3,'Value', 0);
    set(handles.tg4,'Value', 0);
    
    errordlg(['Video dir path: ', mocapstruct.cameradirectory{2}, ' is not valid. Choose a path in Setting Menu.'], ...
        'Video Directory Path Check');
    
    return;
else
    videofolder = mocapstruct.cameradirectory{2};
end

matchedindex = mocapstruct.matched_frames_aligned{2}(frame_inds);

h1 = handles.a1;
h2 = handles.a2;

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
set(h1,'Color','k');

xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1));
yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2));
zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3));
line(h1,xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);

if frame_inds ==1

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

oldframe = 0;
%base_time = datenum(mocapstruct.mocapfiletimes{1});
base_time = datenum(0);

for lk = reshape(frame_inds,1,[])

    cla(h1);

    %% Plot the marker movie
    set(handles.t4,'String',num2str(lk));
    ind_to_plot = lk;

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

    %% plot the video movie    
    
    %videoframe = matchedindex(find(frame_inds==lk)); % ## Not Needed
    videoframe = matchedindex;
    
    framenumber = mod(videoframe, 3500);

    if ( videoframe > 0 && oldframe ~= videoframe )        

        videofile = fullfile(videofolder, strcat(num2str(videoframe-framenumber),'.mp4'));
        
        if ~exist(videofile, 'file')
            set(handles.tg1,'Value', 0);
            set(handles.tg2,'Value', 0);
            set(handles.tg3,'Value', 0);
            set(handles.tg4,'Value', 0);
            errordlg(['Video file: ', videofile, ' is not there.'], ...
                'Video File Check');
            return;
        end

        [imagesout] = MP4Reader(videofile, framenumber);
        imagesc(h2, imagesout.cdata);
        set(h2,'xtick',[],'ytick',[]);

    end

    % Removed since it's not used ##
    % M(find(frame_inds == lk)) =  getframe(h2);

end

set(h1,'Nextplot','add');
set(h2,'Nextplot','add');

end
