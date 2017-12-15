function M = AnimateMarker(mocapstruct,frame_inds,hObject, eventdata, handles)

global CurrentFrame;

h=handles.a1;
% frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));

%M{jjj} = movie;
%Mhere = movie;
frame_last = 0;

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
set(h,'Color','k')

xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1));
yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2));
zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3));
line(h,xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);

if CurrentFrame ==1
    
    axis(h,'manual');
    set(h,'Color','k');
    grid(h,'on');
    set(h,'Xcolor',[1 1 1 ]);
    set(h,'Ycolor',[1 1 1]);
    set(h,'Zcolor',[1 1 1]);
    
    zlim(h,[-150 200]);
    xlim(h,[-400 400]);
    ylim(h,[-400 400]);
    set(h,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[]);
    view(h,[-48, 12]);
    
end


for lk = reshape(frame_inds,1,[])%1:10:10000
    
    cla(h);
    
    set(handles.t4,'String',num2str(lk));
    ind_to_plot = lk;
    
    %% Plot markers that are tracked in the frame
    set(h,'Nextplot','ReplaceChildren');
    handles_here = cell(1,numel(mocapstruct.markernames));
    for jj = 1:numel(mocapstruct.markernames)
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
            if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
                
                xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                handles_here{jj} = line(h,xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',8);
                
                hold(h,'on');
                marker_plot(jj) = 1;
            else
                marker_plot(jj) = 0;
            end
        end
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
                line(h,xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',3);
                
                
            end
        end
    end
    
    title(h,strcat('  Frame: ' ,num2str(lk)),'Color','w')
    
    drawnow;
    hold(h,'off');
    
    frame_last = lk;
    
end
