function createmocapfilestruct(ratname)


if strcmp(ratname,'Vicon8')
mocapfilestruct = struct('PreLesion',[],'UniLesion',[]);
mocapfilestruct.mocapdir = 'E:\Bence\Data\Motionanalysis_captures\Vicon8\';
mocapfilestruct.PreLesion.days = {'20170816\Generated_C3D_files\','20170817\Generated_C3D_files\','20170818\Generated_C3D_files\',...
    '20170820\Generated_C3D_files\',...
    '20170821\Generated_C3D_files\','20170822\Generated_C3D_files\','20170823\Generated_C3D_files\',... day7
    };
mocapfilestruct.UniLesion.days = {'20170824\Generated_C3D_files\','20170825\Generated_C3D_files\',...%8,9
    '20170827\Generated_C3D_files\','20170828\Generated_C3D_files\',...%10,11
    '20170830\Generated_C3D_files\','20170831\Generated_C3D_files\','20170901\Generated_C3D_files\','20170904\Generated_C3D_files\'};

%% specific markersets for each day
    [markercolor,links] = day_loading_header('vicon8_20marker');

filestruct_conds = fieldnames(mocapfilestruct);
for ll = 1:numel(filestruct_conds)
    if isfield(mocapfilestruct.(filestruct_conds{ll}),'days')
       days_here = mocapfilestruct.(filestruct_conds{ll}).days;
   mocapfilestruct.(filestruct_conds{ll}).links = cell(1,numel(days_here));
      mocapfilestruct.(filestruct_conds{ll}).markercolor = cell(1,numel(days_here));
for jj = 1:numel(days_here)
     mocapfilestruct.(filestruct_conds{ll}).links{jj} = links;
      mocapfilestruct.(filestruct_conds{ll}).markercolor{jj} = markercolor;
end
    end
end
end


if strcmp(ratname,'JDM25')
mocapfilestruct = struct('PreLesion',[],'UniLesion',[]);
mocapfilestruct.mocapdir = 'E:\Bence\Data\Motionanalysis_captures\JDM25\';
mocapfilestruct.PreLesion.days = {'20170913\Generated_C3D_files\','20170914\Generated_C3D_files\','20170915\Generated_C3D_files\',...
    '20170916\Generated_C3D_files\','20170917\Generated_C3D_files\','20170918\Generated_C3D_files\','20170919\Generated_C3D_files\','20170920\Generated_C3D_files\'};


%% colos and links
   [markercolor,links] = day_loading_header('vicon8_20marker');

filestruct_conds = fieldnames(mocapfilestruct);
for ll = 1:numel(filestruct_conds)
    if isfield(mocapfilestruct.(filestruct_conds{ll}),'days')
       days_here = mocapfilestruct.(filestruct_conds{ll}).days;
   mocapfilestruct.(filestruct_conds{ll}).links = cell(1,numel(days_here));
      mocapfilestruct.(filestruct_conds{ll}).markercolor = cell(1,numel(days_here));
for jj = 1:numel(days_here)
     mocapfilestruct.(filestruct_conds{ll}).links{jj} = links;
      mocapfilestruct.(filestruct_conds{ll}).markercolor{jj} = markercolor;
end
end
end
end

%% Vicon3 properties (more complicated)
if strcmp(ratname,'Vicon3')
    
       mocapdir = 'E:\Bence\Data\Motionanalysis_captures\Vicon3\';

mocapfilestruct = struct('early',[],'OneHand',[],'eighteenmarker',[],'eighteenmarkernoheadcap',[],'seventeen',[]);
mocapfilestruct.mocapdir = 'E:\Bence\Data\Motionanalysis_captures\Vicon3\';
mocapfilestruct.early.days = {'20170630\Generated_C3D_files\','20170701\Generated_C3D_files\','20170703\Generated_C3D_files\'};
mocapfilestruct.OneHand.days = {'20170705\Generated_C3D_files\','20170706\Generated_C3D_files\',...
    '20170707\Generated_C3D_files\','20170710\Generated_C3D_files\','20170711\Generated_C3D_files\','20170712\Generated_C3D_files\','20170714\Generated_C3D_files\'};
mocapfilestruct.eighteenmarker.days = {'20170720\Generated_C3D_files\','20170721\Generated_C3D_files\','20170722\Generated_C3D_files\'};
mocapfilestruct.eighteenmarkernoheadcap.days = {'20170725\Generated_C3D_files\','20170726\Generated_C3D_files\','20170727\Generated_C3D_files\'};
mocapfilestruct.seventeen.days = {'20170728\Generated_C3D_files\'};

end

filestruct_conds = fieldnames(mocapfilestruct);
for ll = 1:numel(filestruct_conds)
    if isfield(mocapfilestruct.(filestruct_conds{ll}),'days')
   days_here = mocapfilestruct.(filestruct_conds{ll}).days;
   
   mocapfilestruct.(filestruct_conds{ll}).mocapfiles = cell(1,numel(days_here));
      mocapfilestruct.(filestruct_conds{ll}).mocapdatecreate = cell(1,numel(days_here));
      mocapfilestruct.(filestruct_conds{ll}).mocapserialtimes = cell(1,numel(days_here));
      mocapfilestruct.(filestruct_conds{ll}).threshcrossings = cell(1,numel(days_here));
      mocapfilestruct.(filestruct_conds{ll}).numframes = cell(1,numel(days_here));
      mocapfilestruct.(filestruct_conds{ll}).missingtimes = cell(1,numel(days_here));

    for jj = 1:numel(days_here)
        
       %   filenames_here1 = dir(strcat(mocapfilestruct.mocapdir,days_here{jj},'*nolj*recording*.c3d'));
           filenames_here = dir(strcat(mocapfilestruct.mocapdir,days_here{jj},'*recording*nolj*.c3d'));
       %filenames_here = cat(1,filenames_here1,filenames_here2);
       
       %% get unique filenames
       tags_unique = cell(1,numel(filenames_here));
              script_end_unique = cell(1,numel(filenames_here));

           if numel(filenames_here)
for kk = 1:numel(filenames_here)
       tags_unique{kk} = strrep(strrep(strrep(strrep(strrep(filenames_here(kk).name,'nolj',''),...
           'Recording',''),'_',''),'.c3d',''),'day','');
       index = find(isletter(tags_unique{kk}), 1);
              index2 = find(isletter(tags_unique{kk}),1, 'last');

script_end_unique{kk}  = tags_unique{kk}(index:index2);
end
           end
           
           %save the unique endings
           mocapfilestruct.(filestruct_conds{ll}).day_conds{jj} = unique(script_end_unique);
           
           [c,ia,ic] = unique(tags_unique,'legacy');
           unique_filenames_here = cell(1,numel(c));
           
           for ihere = 1:numel(c)
               filename_inds = find(ic == ihere);
               size_fname = zeros(1,numel(filename_inds));
               for lk = 1:numel(filename_inds)
                   size_fname(lk) = numel(filenames_here(filename_inds(lk)).name);
               end
               [~,indmax] = max(size_fname);
               unique_filenames_here{ihere} = filenames_here(filename_inds(indmax)).name;
               
           end
           %% only loop over the unique files
        fprintf('on day %s number of files %f \n',days_here{jj},numel(unique_filenames_here));
         filepath_array = [];
         filepath_ind = 1;
if numel(unique_filenames_here)

    for kk = 1:numel(unique_filenames_here)
        if (numel(strfind(unique_filenames_here{kk},'nolj')))
            filepath_array{filepath_ind} = strcat(mocapfilestruct.mocapdir,days_here{jj},unique_filenames_here{kk});
            filepath_ind = filepath_ind+1;
        end
    end
    
    [mocapfilestruct.(filestruct_conds{ll}).mocapdatecreate{jj},...
        mocapfilestruct.(filestruct_conds{ll}).mocapserialtimes{jj},...
        mocapfilestruct.(filestruct_conds{ll}).mocapfiles{jj}] = sort_mocap_files(filepath_array,mocapfilestruct.mocapdir);
    
    [~,mocapfilestruct.(filestruct_conds{ll}).threshcrossings{jj},...
        mocapfilestruct.(filestruct_conds{ll}).numframes{jj},...
        mocapfilestruct.(filestruct_conds{ll}).missingtimes{jj}] = get_analogactive_files(mocapfilestruct.(filestruct_conds{ll}).mocapfiles{jj});
end
    end
    end
    
    %% get the file times and sort each
     %   [mocap_datecreate,sorted_mocap_serialtimes,filepath_array_sorted ] = sort_mocap_files(filepath_array,mocapdir);

    
end


save(strcat(mocapfilestruct.mocapdir,'mocapfilestruct_',ratname,'_.mat'),'mocapfilestruct');


for ll = 1:numel(filestruct_conds)
     if isfield(mocapfilestruct.(filestruct_conds{ll}),'days')
         fprintf('For type %s nhours %f nframes %f \n',(filestruct_conds{ll}),...
             sum(cellfun(@sum,mocapfilestruct.(filestruct_conds{ll}).numframes))./(300*3600),...
             sum(cellfun(@sum,mocapfilestruct.(filestruct_conds{ll}).numframes)));
     end
end
%Vicon3 ~ 150 hrs, 162 million frames


