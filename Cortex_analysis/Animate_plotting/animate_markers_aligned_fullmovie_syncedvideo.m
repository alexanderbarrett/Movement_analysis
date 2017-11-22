function M = animate_markers_aligned_fullmovie_syncedvideo(mocapstruct,frame_inds,videofolder,matchedindex)
%matlab_fr = 10;
h=figure(370)
% frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));

%M{jjj} = movie;
%Mhere = movie;
frame_last = 0;

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
subplot(2,1,1)
    set(h,'Color','k')
    plot3( squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1)),...
        squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2)),...
        squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3)),'o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6)
    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
    grid on;
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    
    zlim([-100 150])
    xlim([-300 300])
    ylim([-300 300])
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
        view([-22, 12]);
ax_h = gca;





oldvidobj = [];
oldfile = 'oldfile.mp4';
oldframe = 0;
fighandle = [];
        base_time = datenum(mocapstruct.mocapfiletimes{1});

for lk = reshape(frame_inds,1,[])%1:10:10000
    

    %% Plot the marker movie
    ind_to_plot = lk;
subplot(2,1,1)
     cla;

%% Plot markers that are tracked in the frame
  set(ax_h,'Nextplot','ReplaceChildren');
    handles_here = cell(1,numel(mocapstruct.markernames));
    for jj = 1:numel(mocapstruct.markernames)
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
        if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
%             handles_here{jj} = plot3( squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1)),...
%                 squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2)),...
%                 squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3)),'o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',8);
%             hold on
  
            xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                handles_here{jj} = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',8);
                %hold on
                
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
%                 plot3( [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
%                     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ],...
%                     [ squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2))...
%                     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))],...
%                     [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3))...
%                     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))],'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'Linewidth',3);
%         
             
                  xx = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
                yy = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
                zz = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
                line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',3);
                
                
                
            end
        end
    end

        title(strcat('  Frame: ' ,datestr(base_time+lk./(mocapstruct.fps*60*60*24))),'Color','w')
    
    
    drawnow
    hold off
  
    frame_last = lk;
    %clf
   
    
    %% plot the video movie
    videoframe = matchedindex(find(frame_inds==lk));
    framenumber = mod(videoframe,3500);
    if (videoframe>0 && oldframe ~=videoframe)
    videofile = strcat(videofolder,num2str(videoframe-framenumber),'.mp4');
    
       subplot(2,1,2)
    [imagesout,oldvidobj,oldfile,fighandle] = mpfour_reader_singleframe(videofile,framenumber,oldvidobj,oldfile,oldframe,fighandle);
 imagesc(imagesout(1).cdata)
    oldframe = framenumber;
    set(gca,'xtick',[],'ytick',[])
%axis off
    end
    
    M(find(frame_inds == lk)) =  getframe(gcf);
end
set(ax_h,'Nextplot','add');

end