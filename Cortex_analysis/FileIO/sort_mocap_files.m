function [mocap_datecreate,sorted_mocap_serialtimes,filepath_array_sorted] = sort_mocap_files(filepath_array,mocapdir)


%% get the date/time of the capture files
%filepaths_capfiles = rdir(strcat(mocapdir,'*',filesep,'*',num2str(EFiles(f)),'*','.anb'));
if (~strcmp(class(filepath_array),'cell'))
    filepath_array_temp = cell(1,1);
    filepath_array_temp{1} = filepath_array;
    filepath_array = filepath_array_temp;
end
filepath_use_1 = strrep(strrep(filepath_array{1},'Generated_C3D_files\',''),'_nolj.c3d','.anb');

regans =  regexp(filepath_use_1,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match');

datetime_array_cap = cell(1,numel(filepath_array));
serialtime_array_cap = zeros(1,numel(filepath_array));
fileending_array_cap = zeros(1,numel(filepath_array));


% change to the current directory to run the dir commany
cd(strcat(mocapdir,regans{2}));
for kk =1:numel(filepath_array)
    if numel(strfind(filepath_array{kk},'Generated_C3D_files\nolj_'))
    filepath_use = strrep(strrep(filepath_array{kk},'Generated_C3D_files\nolj_',''),'_nolj.c3d','.anb');
    else
            filepath_use = strrep(strrep(filepath_array{kk},'Generated_C3D_files\',''),'_nolj.c3d','.anb');
    end
    
    stringcmd = char(strcat('dir',{' '},filepath_use ));
    [dum,str] = dos(stringcmd);
    c = textscan(str,'%s');
    t = (datetime(strcat(c{1}{15},'/',c{1}{16},':00/',c{1}{17}),'InputFormat','MM/dd/yyyy/hh:mm:ss/a'));
    stringcmd;
    %  t
    %t = double(t);
    %datehash = t(3)+t(4)./24+t(5)./(24*60);
    
    datetime_array_cap{kk} = t;
    serialtime_array_cap(kk) = datenum(t);
    
    %% also search to get the fileending and make sure this agrees with the fileending
    regans =  regexp(filepath_use ,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match');
    

end

v = cellfun(@datestr,datetime_array_cap,'UniformOutput', false);

[aa,bb] = sort(serialtime_array_cap,'ASCEND');
mocap_datecreate = v(bb);
sorted_mocap_serialtimes= serialtime_array_cap(bb);
filepath_array_sorted = filepath_array(bb);
