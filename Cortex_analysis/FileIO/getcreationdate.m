function [datetime_array,serialtime_array] = getcreationdate(filepath_use,mocapdir,dayuse)
 
cd(strcat(mocapdir,dayuse));

stringcmd = char(strcat('dir',{' '},filepath_use ));
    [dum,str] = dos(stringcmd);
    c = textscan(str,'%s');
    t = (datetime(strcat(c{1}{15},'/',c{1}{16},':00/',c{1}{17}),'InputFormat','MM/dd/yyyy/hh:mm:ss/a'));
    stringcmd;
    %  t
    %t = double(t);
    %datehash = t(3)+t(4)./24+t(5)./(24*60);
    
    datetime_array = t;
    serialtime_array = datenum(t);
    
    %% also search to get the fileending and make sure this agrees with the fileending
    regans =  regexp(filepath_use ,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match');
end