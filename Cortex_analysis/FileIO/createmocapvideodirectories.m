function createmocapvideodirectories(ratname)


if strcmp(ratname,'Vicon8')
mocapfilestruct = struct('PreLesion',[],'UniLesion',[]);
mocapfilestruct.mocapdir = '\\140.247.178.37\Motionanalysis_captures\Vicon8\';
mocapfilestruct.PreLesion.days = {'20170816\','20170817\','20170818\',...
    '20170820\',...
    '20170821\','20170822\','20170823\',... day7
    };
mocapfilestruct.UniLesion.days = {'20170824\','20170825\',...%8,9
    '20170827\','20170828\',...%10,11
    '20170830\','20170831\','20170901\','20170904\'};
end


if strcmp(ratname,'JDM25')
mocapfilestruct = struct('PreLesion',[],'UniLesion',[]);
mocapfilestruct.mocapdir = 'Y:\Jesse\Data\Motionanalysis_captures\JDM25\';
mocapfilestruct.PreLesion.days = {'20170913\','20170914\','20170915\',...
    '20170916\','20170917\','20170918\','20170919\','20170920\','20170921\','20170922\','20170923\'};
mocapfilestruct.UniLesion.days = {'20170924\','20170925\','20170926\','20170927\','20170928\','20170929\','20171002\','20171003\'};
mocapfilestruct.BiLesion.days = {'20171004\','20171005\','20171006\','20171007\','20171008\','20171009\','20171010\','20171011\','20171012\','20171013\'};

end


if strcmp(ratname,'JDM32')
mocapfilestruct = struct('PreLesion',[],'UniLesion',[],'BiLesion',[]);
mocapfilestruct.mocapdir = 'Y:\Jesse\Data\Motionanalysis_captures\JDM32\';
mocapfilestruct.PreLesion.days = {'20171020\','20171022\','20171023\',...
    '20171024\','20171025\','20171026\','20171027\','20171027_2\','20171028\','20171029\','20171030\'};
mocapfilestruct.UniLesion.days = {'20171031\','20171031_2\','20171101\','20171101_2\','20171102\','20171103\',...
    '20171103_2\','20171104\','20171105\','20171106\','20171107\'};
mocapfilestruct.BiLesion.days = {'20171110\','20171113\','20171113_2\','20171114\','20171115\','20171116\','20171117\','20171118\'};

end


if strcmp(ratname,'Vicon3')
    
  
    mocapdir = ' \\140.247.178.37\Motionanalysis_captures\Vicon3\';

mocapfilestruct = struct('early',[],'OneHand',[],'eighteenmarker',[],'eighteenmarkernoheadcap',[],'seventeen',[]);
mocapfilestruct.mocapdir = 'E:\Bence\Data\Motionanalysis_captures\Vicon3\';
mocapfilestruct.early.days = {'20170630\','20170701\','20170703\'};
mocapfilestruct.OneHand.days = {'20170705\','20170706\','20170707\','20170710\','20170711\','20170712\','20170714\'};
mocapfilestruct.eighteenmarker.days = {'20170720\','20170721\','20170722\'};
mocapfilestruct.eighteenmarkernoheadcap.days = {'20170725\','20170726\','20170727\'};
mocapfilestruct.seventeen.days = {'20170728\'};

end


filestruct_conds = fieldnames(mocapfilestruct);

for ll = 1:numel(filestruct_conds)
    if isfield(mocapfilestruct.(filestruct_conds{ll}),'days')
   days_here = mocapfilestruct.(filestruct_conds{ll}).days;

    for jj = 1:numel(days_here)

        
       %   filenames_here1 = dir(strcat(mocapfilestruct.mocapdir,days_here{jj},'*nolj*recording*.c3d'));
           filenames_here = dir(strcat(mocapfilestruct.mocapdir,days_here{jj},'*recording*.cap'));
           
           
       %filenames_here = cat(1,filenames_here1,filenames_here2);
      % filepath_array{1} = filenames_here(2).name;
    


    
    
       %% get unique filenames
       tags_unique = cell(1,numel(filenames_here));
              script_end_unique = cell(1,numel(filenames_here));

           if numel(filenames_here)
for kk = 1:numel(filenames_here)
       tags_unique{kk} = strrep(filenames_here(kk).name,'.cap','');
      % tags_unique{kk} = filenames_here(kk).name;
       index = find(isletter(tags_unique{kk}), 1);
              index2 = find(isletter(tags_unique{kk}),1, 'last');
if (numel(strfind(tags_unique{kk},'636'))==0)
script_end_unique{kk}  = tags_unique{kk}(index:index2);
else
    %if doing ephys keep the number
    end_no = min(index2+10,numel(tags_unique{kk}));
    script_end_unique{kk}  = tags_unique{kk}(index:( end_no));
end
end
           end
           
           %save the unique endings
           unique_folder_names = unique(script_end_unique);
           
           %% move the folders and transfer the files
           videosubscripts = {'CameraL','CameraU','CameraR'};
               videopathbase =  '\\MOCAPRIG\MOCAP_training\Testing\';

               %get the creation date of a file with the unique name

               fprintf('Copying Video directories and making folders \n')
           for kjk = 1:numel(unique_folder_names)
                         uniquefiles = find(cellfun(@numel,strfind(tags_unique,unique_folder_names{kjk}))==1);
                         serialarray = zeros(1,numel(uniquefiles));
                         serialarray = double(serialarray);
                         for lkk = 1:numel(uniquefiles)
                             % binary file isn't modified
                    [there,serialarray(lkk)] =  getcreationdate(strcat(tags_unique{uniquefiles(lkk)},'.anb'),mocapfilestruct.mocapdir,days_here{jj});
                         end
               serialt = min(serialarray);
               
              mkdir(strcat(mocapfilestruct.mocapdir,days_here{jj},unique_folder_names{kjk} ));
              for mm =1:numel(videosubscripts)
                    outputdir = strcat(mocapfilestruct.mocapdir,days_here{jj},unique_folder_names{kjk},filesep,videosubscripts{mm}) ;

                            mkdir(outputdir);
                            videopathFolder = strcat(videopathbase,videosubscripts{mm});
                            
                                 [videofolders,nameFolds] = getfoldernumbers(videopathFolder);
                                 
                                 %% there is a 3 hr timezone offset and a 1-2 min grace period
    timedifferences =  (-10^7*3600*2.9+double(convertMATLABtime((serialt))))-videofolders;
    folderind = find(timedifferences>0,1,'last');
    
    %there is a 3 hr timezone offset
    
    
    %% if less than a 2 hr gap
   if (timedifferences(folderind)< 10^9*3600*2)
       
  videosourcehere = strcat(videopathFolder,filesep,nameFolds(folderind));
  outputfolder = strcat(outputdir,filesep,nameFolds(folderind));
  if (~exist(outputfolder{1},'dir'))
      fprintf('Making output directory for day %f condition %f \n',jj,ll)
  mkdir(outputfolder{1});
  copyfile(videosourcehere{1},outputfolder{1})
  end
              end

              end
           end
        
           
           %mocap videos
           %
           
    end
    end
    
    %% get the file times and sort each
     %   [mocap_datecreate,sorted_mocap_serialtimes,filepath_array_sorted ] = sort_mocap_files(filepath_array,mocapdir);

    
end
