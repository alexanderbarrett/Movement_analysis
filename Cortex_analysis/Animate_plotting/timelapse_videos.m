mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots_videoexamples',filesep);
mkdir(savedirectory);

close all;
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);

%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server

movie_matrix = cell(2,max(numel(mocapfilestruct.PreLesion.days),numel(mocapfilestruct.UniLesion.days)));

for taghere = 1:2
    if (taghere == 1)
descriptor_struct_1.cond = 'PreLesion';
    else
      descriptor_struct_1.cond = 'UniLesion';
  
    end
descriptor_struct_1.tag = 'overnight';

for kk=1:numel(mocapfilestruct.(descriptor_struct_1.cond).days)
       descriptor_struct_1.day = kk;

if sum(ismember(mocapfilestruct.(descriptor_struct_1.cond).day_conds{kk},'overnight'))>0
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
good_inds_2 = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},'morning')));

good_inds = strcat(good_inds,good_inds_2);
good_inds = good_inds(1:4);
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
 mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);
    
%[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_videotest',mocapmasterdirectory);

[mocapstruct] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes);

[modular_cluster_properties] = get_modularclustproperties(mocapstruct);

v = VideoWriter(strcat(savedirectory,'sbys_movie_',num2str(taghere),'lesion',num2str(kk)),'MPEG-4');
                videolength = 100*10*60;
                open(v)                
                
                time_subset = intersect(mocapstruct.move_frames_fast,modular_cluster_properties.clustering_inds_agg{8});
                
                
                time_subset_contig = find_contig_blocks(time_subset,300);
                
                figure(333)
                plot(medfilt1(mocapstruct.markers_preproc.SpineL(time_subset,3),6));
                
            %    inds = intersect(modular_cluster_properties_social.clustering_inds_agg{2},,
movie_matrix{taghere ,kk} = animate_markers_aligned_fullmovie(mocapstruct,time_subset((1:30:min(numel(time_subset),videolength))));
   writeVideo(v,movie_matrix{taghere ,kk})
                                    close(v)
end
end
end



plot_simul = 1;
if (plot_simul)
    for taghere = 1

    nrows = 2;
    ncols = 3;
    
    %  cluster_numbers = [2,42,51,55,59,124,73,69,65,62,183,113,112,164,158,175,15,183];
    %             cluster_numbers = [6,50,11,43,35,33,29,22,16,13,12,10];
    cluster_numbers = [2,3,4,5,6,7];
    
    %      good_clusters(1:16);
    %            cluster_names = {'tap/drop','groom hands','odor sample','up/down/up',....
    %                'under and up','tap and rise','sniff and explore','left tap','lick',...
    %                'tap/drop','hi sample','tap','lick and up','med sample','lever sample','low sample'};
    %                     cluster_names = {'hi sample','sniff/explore','mid sniff','mid sample','low sample','head crane','reach hisample','low sniff',...
    %                         'eating','sniff low','down and up','sniff and shake','sniff and explore'};
    cluster_names = {'Day 2','Day 3','Day 4','Day 5','Day 6', 'Day 7'};
    
    
    %num2str(cluster_numbers');
    fighand =   figure(388);
    set(fighand,'Color','k')
    set(fighand,'Position',[100 100 1100 1100])
    for lk = 1:1000 %numel(movie_output{jjj})
        for ll = 1:nrows*ncols  %numel(cluster_numbers)
            mov_ind = cluster_numbers(ll);
            subplot_tight(nrows,ncols, ll)
            movie_size = numel(movie_matrix{mov_ind});
            frame_use = mod(lk,movie_size);
            if (frame_use==0)
                frame_use = 1;
            end
            
         %   frameuse = mod(
            
            imshow(movie_matrix{taghere,mov_ind}(frame_use).cdata,movie_matrix{taghere,mov_ind}(frame_use).colormap)
            title(cluster_names{ll},'Color','w')
        end
        M_full(lk) =      getframe(gcf);
        
    end
    
    v = VideoWriter(strcat(savedirectory,'aggregate_movie_lesion',num2str(taghere)),'MPEG-4');
    open(v)
    writeVideo(v, M_full)
    close(v)
    clear M_full
    end
end
                

%% also save an accelerated video
